"""oalib.py
Library for online analysis in the Baccus Lab.
(C) 2017 Benjamin Naecker bnaecker@stanford.edu
"""

import abc
import warnings

import numpy as np
import h5py
import matplotlib.pyplot as plt
from scipy import signal
import pyret

class OnlineAnalysis(metaclass=abc.ABCMeta):
    '''Abstract base class implementing the interface for online 
    analysis methods.
    '''

    def __init__(self, stimfile, frame_rate, length, channels):
        '''Construct an OnlineAnalysis object.

        This abstract base class defines the interface for all online analysis
        objects, and implements a few functions that help in that analysis.

        Parameters
        ----------

        stimfile : str
            A filename of the HDF5 file containing the stimulus. This file
            must have a dataset called '/stimulus', which should be shaped
            as `(t, x, y)`, i.e., time by space by space. Either or both
            of the spatial dimensions may be omitted, for example for a simple
            full-field temporal stimulus.

        frame_rate : float
            The effective frame rate of the stimulus, in Hz. If the value of
            0.0 or None is passed, then the `/stimulus` dataset must contain an 
            attribute called `frame-rate`, which will be read and used as the 
            effective rate.

        length : float
            The length of the analysis, in seconds. This is usually something
            like the length of an online receptive field.

        channels : array_like of int
            A list of channels on which the analysis will be performed. This
            is usually a single channel, but for more exotic analyses, it may
            contain multiple channels.

        Helpers
        -------

        There are several useful helper functions defined in the base class. This
        allows subclasses to easily call them to aid analysis.

            - `frame_index_from_time(t)`
                Given a floating-point time value, return the stimulus frame 
                number on the screen at that time.

            - `frame_from_time(t)`
                Given a floating-point time value, return the stimulus frame 
                itself on the screen at that time.

            - `frame_times(start, stop)`
                Return an array giving the time of each stimulus frame within
                the given start and stop times.

            - `frames_from_time(start, stop)`
                Return the full stimulus frames between the start and stop
                times given.

        There are also several public data attributes which can be used to aid
        in analysis:

            - `stimulus`
                The HDF5 stimulus dataset.

            - `frame_rate`
                The effective frame rate of the stimulus.

            - `figure`
                The matplotlib figure into which the analysis will be plotted.

            - `current_result`
                The current result of the online analysis.

            - `current_time`
                The time of the most recent data received.

            - `length`
                The length in seconds of the analysis. For example, the length
                of an online receptive field estimate.

            - `channels`
                The channels in the data to which the analysis is applied.

        Subclassing
        -----------

        The main (abstract) methods that must be implemented in subclasses are:

            - `update()`
                Update the current analysis result from a new chunk of data.

            - `plot()`
                Plots the current result of the analysis. This need not actually
                do anything.

            - `finalize()`
                Called after the end of analysis to finalize both the result
                itself and the plot of it.

        The current result of the analysis should be stored in the `current_result`
        data attribute. Subclass constructors should fill this field with
        whatever value is needed. For example, for a subclass computing a receptive
        field, this would be initialized with a zeros array of the right size.

        See the docstrings for those methods for details on how to re-implement
        them in subclasses.
        '''
        self._stimfile_name = stimfile
        self.frame_rate = frame_rate
        self.length = length
        self.channels = channels

        self._init_figure()
        self._open_stimfile()
        self.current_result = None
        self.current_time = 0.0


    def __str__(self):
        return 'OnlineAnalysis(stimfile={}, length={:0.2f}, channels={})'.format(
                self._stimfile_name, self.length, self.channels)


    def __del__(self):
        if self.figure is not None:
            plt.close(self.figure)
        try:
            self._stimfile.close()
        except:
            pass


    @abc.abstractmethod
    def update(self, dataframe):
        '''Updates the current result of the analysis with a new frame of data.

        Parameters
        ----------

        dataframe : bldsclient.DataFrame
            A new frame of data in the recording. This is defined in the `bldsclient`
            module, and is just an ndarray containing the data, plus a `start` and
            `stop` time for the chunk of data.

        This method must update the `current_result` data attribute of the subclass,
        using whatever analyses are appropriate for that subclass. For example, for
        an online reverse-correlation, this would compute the reverse-correlation 
        between the new data and the corresponding stimulus, and then just add this
        to the `current_result` field.
        '''
        pass


    @abc.abstractmethod
    def plot(self):
        '''Plot the current result of the analysis.

        This method should be reimplemented in subclasses and plot whatever feature
        of the current result is appropriate. For example, for an online temporal
        receptive field analysis, this would just plot that RF.

        Subclass implementations should plot data into the figure defined by the
        `figure` data attribute. If that figure is None, it means that it was 
        closed, probably by the user at some point. The subclass implementations
        must check if the figure exists, and recreate it and any children before
        doing any actual plotting.
        '''
        pass


    @abc.abstractmethod
    def finalize(self):
        '''Finalize the result of the analysis.

        This method should be reimplemented after an online analysis has completed,
        and should do whatever computation and plotting is required to finish the
        analysis. For example, an online receptive field estimate might be
        smoothed before being replotted.
        '''
        pass


    def _open_stimfile(self):
        '''Open the stimulus file and check it has a frame rate.'''
        self._stimfile = h5py.File(self._stimfile_name, 'r')
        if 'stimulus' not in self._stimfile:
            self._stimfile.close()
            raise ValueError('''The requested stimulus file '{}' does
                    not contain a '/stimulus' attribute.'''.format(
                    self._stimfile_name))
        self.stimulus = self._stimfile['stimulus']

        if self.frame_rate is None or self.frame_rate == 0.0:
            if 'frame-rate' not in list(self._stimfile['stimulus'].attrs):
                self._stimfile.close()
                raise ValueError('''No frame rate was specified, but the
                        requested stimulus file '{}' does not contain one.
                        '''.format(self._stimfile_name))
            self.frame_rate = self.stimulus.attrs.get('frame-rate')

        self._stim_time = np.arange(self.stimulus.shape[0]) / self.frame_rate


    def _init_figure(self):
        self.figure = plt.figure()
        self._cid = self.figure.canvas.mpl_connect('close_event', self._on_figure_close)

    
    def _on_figure_close(self, event):
        '''Small handler function to set the self.figure attribute to None when
        the figure is closed. This is used so that subclasses can learn when 
        the figure has disappeared and they need to re-create it in the `plot()`
        method.
        '''
        self.figure.canvas.mpl_disconnect(self._cid)
        self.figure = None
        self._cid = None


    def frame_index_from_time(self, t):
        '''Return the index of the frame at the given time.

        Parameters
        ----------

        t : float
            The time at which the frame number should be returned.

        Raises
        ------

        If the time is invalid for the stimulus, e.g., beyond the length
        of the stimulus, a ValueError is raised.
        '''
        ix = np.where(self._stim_time >=t)[0]
        if len(ix):
            return ix[0]
        raise ValueError('''The time {} is invalid for the stimulus.
                '''.format(t))


    def frame_from_time(self, t):
        '''Return the stimulus frame at the given time.

        Parameters
        ----------

        t : float
            The time of the frame to be returned.

        Raises
        ------

        If the time is invalid for the stimulus, e.g., beyond the length
        of the stimulus, a ValueError is raised.
        '''
        return self.stimulus[self.frame_index_from_time(t), ...]


    def frame_times(self, start, stop):
        '''Return the time of the stimulus frames between the given values.

        Parameters
        ----------

        start, stop : float
            The start and stop times between which the time of each
            stimulus frame is to be returned.

        Raises
        ------

        If one or more of the given times is invalid for the stimulus, e.g.,
        beyond the length of the stimulus, a ValueError is raised.
        '''
        return self._stim_time[self.frame_index_from_time(start) : 
                self.frame_index_from_time(stop)]


    def frames_from_time(self, start, stop):
        '''Return the stimulus frames between the given times.

        Parameters
        ----------

        start, stop : float
            The start and stop times between which the stimulus frames
            are returned.

        Raises
        ------

        If one or more of the given times is invalid for the stimulus, e.g.,
        beyond the length of the stimulus, a ValueError is raised.
        '''
        return self.stimulus[self.frame_index_from_time(start) : 
                self.frame_index_from_time(stop), ...]

class OnlineReceptiveField(OnlineAnalysis):
    '''Abstract base class for computing receptive fields online.'''

    def __init__(self, stimfile, frame_rate, length, channel, 
            filter_sample_rate=1000.0):
        '''Construct an OnlineReceptiveField class.

        Parameters
        ----------
        
        stimfile : str
            Filename of an HDF5 file containing the stimulus. 
            Refer to the OnlineAnalysis base class for more details.

        frame_rate : float
            The effective frame rate of the stimulus.
            Refer to the OnlineAnalysis base class for more details.

        length : float
            Length of the receptive field, in seconds.
            Refer to the OnlineAnalysis base class for more details.

        channel : int
            An integer or array_like with one integer element, giving
            the channel to use for computing the filter. This is the
            only argument which differs from the base class implementation.

        filter_sample_rate : float, optional
            The desired sample rate of the resulting receptive field, in Hz.
            Defaults to 1000.0. The data and stimulus will both be resampled
            to this rate.

        Notes
        -----

        This class provides the infrastructure for plotting an RF
        online. It does not actually compute an RF, which is implemented
        in the subclasses, e.g., OnlineReverseCorrelation.
        '''
        # Verify input is a single channel
        try:
            nchannels = len(channel)
            channel = int(channel[0])
        except TypeError:
            channel = int(channel)
        super().__init__(stimfile, frame_rate, length, channel)

        # Initialize receptive field
        self.filter_sample_rate = filter_sample_rate
        self.filter_length = int(np.round(self.length * self.filter_sample_rate))
        self.current_result = np.zeros((self.filter_length,) + 
                self.stimulus.shape[1:])

        self._init_figure_window()


    def __str__(self):
        return 'OnlineReceptiveField(stimfile={}, length={:0.2f}, channels={})'.format(
                self._stimfile_name, self.length, self.channels)


    def plot(self):
        '''Plot 1, 2, or 3D receptive field with the current result.

        This method will plot either a single temporal receptive field,
        a 1D spatiotemporal receptive field, or a decomposed 2D spatiotemporal
        receptive field as a temporal and spatial kernel.
        '''

        # Recreate figure if needed
        if self.figure is None:
            self._init_figure()
            self._init_figure_window()

        # Plot 1D (pure temporal RF)
        if self.stimulus.ndim == 1:
            self._temporal_rf_line.set_ydata(self.current_result)
            mn = self.current_result.min()
            mn -= 0.1 * abs(mn)
            mx = self.current_result.max()
            mx += 0.1 * abs(mx)
            self._axes[0].set_ylim(mn, mx)

        # Plot 2D spatiotemporal RF
        elif self.stimulus.ndim == 2:
            self._spatiotemporal_rf_image.set_data(self.current_result.T)
            self._spatiotemporal_rf_image.set_clim(
                    self.current_result.min(),
                    self.current_result.max())

        # Plot decomposed spatial and temporal kernels
        else:
            self._temporal_rf_line.set_ydata(self._temporal_rf)
            mn = self._temporal_rf.min()
            mn -= 0.1 * abs(mn)
            mx = self._temporal_rf.max()
            mx += 0.1 * abs(mx)
            self._axes[0].set_ylim(mn, mx)
            self._spatial_rf_image.set_data(self._spatial_rf)
            self._spatial_rf_image.set_clim(self._spatial_rf.min(),
                    self._spatial_rf.max())

        # Update time 
        self.figure.suptitle(self.current_title)

        # Force update of the plot
        plt.pause(0.01)


    def _init_figure_window(self):
        '''Set up a figure with axes and lines appropriate for
        the dimensionality of the computed RF.
        '''
        self._axes = []

        # Single axis and line for the temporal RF
        if self.stimulus.ndim == 1:
            self._axes.append(self.figure.add_subplot(111))
            self._axes[0].set_title('Temporal RF')
            self._axes[0].set_xlabel('Time (s)')
            self._axes[0].set_yticks([])
            self._axes[0].set_xticks((0, self.length / 2, 
                    self.length))

            xx = np.linspace(0, self.length, self.filter_length)
            self._temporal_rf_line = self._axes[0].plot(xx,
                    self.current_result, linewidth=3)[0]

        # Single image showing 2D spatiotemporal RF
        elif self.stimulus.ndim == 2:
            self._axes.append(self.figure.add_subplot(111))
            self._axes[0].set_title('Spatiotemporal RF')
            self._axes[0].set_xlabel('Time (s)')
            self._axes[0].set_ylabel('Space')
            self._axes[0].set_xticks((0, 
                    self.filter_length / 2, 
                    self.filter_length))
            self._axes[0].set_xticklabels(('0', 
                    str(self.length / 2), 
                    str(self.length)))

            self._spatiotemporal_rf_image = self._axes[0].imshow(
                    self.current_result.T, origin='lower')

        # Two axes, one each for spatial and temporal kernels
        else:
            self._axes.append(self.figure.add_subplot(121))
            self._axes[0].set_title('Temporal RF')
            self._axes[0].set_xlabel('Time (s)')
            self._axes[0].set_yticks([])
            self._axes[0].set_xticks((0, self.length / 2, 
                    self.length))

            xx = np.linspace(0, self.length, self.filter_length)
            self._temporal_rf_line = self._axes[0].plot(xx,
                    np.zeros((self.filter_length,)), linewidth=3)[0]

            self._axes.append(self.figure.add_subplot(122))
            self._axes[1].set_title('Spatial RF')

            self._spatial_rf_image = self._axes[1].imshow(
                    np.zeros(self.stimulus.shape[1:]), origin='lower')

            # Make dummy kernels
            self._spatial_rf, self._temporal_rf = pyret.filtertools.decompose(
                    self.current_result)

        self.current_title = ''


class OnlineReverseCorrelation(OnlineReceptiveField):
    '''Class performing online reverse-correlation.'''

    def __init__(self, stimfile, frame_rate, length, channel):
        '''Construct an OnlineReverseCorrelation object.

        Parameters
        ----------
        
        stimfile : str
            Filename of an HDF5 file containing the stimulus. 
            Refer to the OnlineAnalysis base class for more details.

        frame_rate : float
            The effective frame rate of the stimulus.
            Refer to the OnlineAnalysis base class for more details.

        length : float
            Length of the receptive field, in seconds.
            Refer to the OnlineAnalysis base class for more details.

        channel : int
            An integer or array_like with one integer element, giving
            the channel to use for computing the filter. This is the
            only argument which differs from the base class implementation.

        Notes
        -----

        The OnlineReverseCorrelation class reimplements the `update()`
        method to compute the reverse-correlation between the stimulus
        and a chunk of the response. This can be used to compute an
        online receptive field estimate when recording intracellularly.
        '''
        super().__init__(stimfile, frame_rate, length, channel)


    def __str__(self):
        return 'OnlineReverseCorrelation(stimfile={}, length={:0.2f}, channels={})'.format(
                self._stimfile_name, self.length, self.channels)


    def update(self, dataframe):
        '''Update the current result with new data.

        Parameters
        ----------
        
        dataframe : DataFrame object
            The chunk of data used to update the online receptive
            field estimate.

        Raises
        ------

        If the time limits of the data frame are outside the bounds of 
        the stimulus (e.g., the recording is longer than the stimulus),
        a ValueError is raised.

        Notes
        -----

        This method implements reverse-correlation between the response
        and the corresponding stimulus chunk. If the stimulus is 3D
        (two spatial dimensions and one temporal), it is decomposed
        into a spatial and temporal kernel as well.
        '''
        self.current_time = dataframe.stop

        # Select channel data
        data = dataframe[self.channels, :].astype(np.double)

        # Compute sampling rate and number of samples per frame
        sample_rate = int(np.round(dataframe.nsamples() / 
                (dataframe.stop - dataframe.start)))
        nsamples_per_frame = int(np.round(sample_rate / 
                self.filter_sample_rate))

        # Detrend data by low-pass filtering and subtracting the
        # filtered signal from the original
        cutoff = 0.5 # Hz
        b, a = signal.butter(3, 0.5 / (sample_rate / 2), btype='low')
        filtered = signal.filtfilt(b, a, data)
        data = (data - filtered)[::nsamples_per_frame]

        # Get frames for this chunk of data and upsample to 1000Hz
        stim = self.frames_from_time(dataframe.start, dataframe.stop)
        factor = int(np.round(self.filter_sample_rate / self.frame_rate))
        stim_upsampled = pyret.stimulustools.upsample(stim, factor)[0]

        # Ensure the stimulus and data have the same size.
        # Repeat the most recent frame at the front of the stimulus.
        if stim_upsampled.shape[0] != data.shape[0]:
            diff = abs(data.shape[0] - stim_upsampled.shape[0])
            stim_upsampled = np.concatenate((stim_upsampled,
                    np.repeat(stim_upsampled[-1, ...].reshape(1, *stim.shape[1:]), 
                    diff, axis=0)), axis=0)
        stim_upsampled -= stim_upsampled.mean(axis=0)

        # Compute reverse correlation, flip along time axis
        self.current_result += pyret.filtertools.revcorr(
                stim_upsampled, data, self.filter_length)[0][::-1, ...]

        # Decompose if necessary
        if self.stimulus.ndim == 3:
            self._spatial_rf, self._temporal_rf = pyret.filtertools.decompose(
                    self.current_result)

        self.current_title = '{:#d}s'.format(int(self.current_time))


    def finalize(self):
        '''Subclass override of the finalize method which does nothing.'''
        pass


class OnlineSpikeTriggeredAverage(OnlineReceptiveField):

    def __init__(self, stimfile, frame_rate, length, channel, spike_method='peaks'):
        '''Construct an OnlineSpikeTriggeredAverage object.

        Parameters
        ----------

        stimfile : str
            Filename of an HDF5 file containing the stimulus. 
            Refer to the OnlineAnalysis base class for more details.

        frame_rate : float
            The effective frame rate of the stimulus.
            Refer to the OnlineAnalysis base class for more details.

        length : float
            Length of the receptive field, in seconds.
            Refer to the OnlineAnalysis base class for more details.

        channel : int
            An integer or array_like with one integer element, giving
            the channel to use for computing the filter. This is the
            only argument which differs from the base class implementation.

        spike_method : str
            The method used to estimate "spike" times from raw voltage data.
            The value of 'peaks' estimates spike times by smoothing the
            voltage data and taking the peak values. A value of 'thresh'
            will estimate spike times by taking upward threshold crossings.
            These should both take about the same amount of time.

        Notes
        -----

        The OnlineSpikeTriggeredAverage class reimplements the `update()`
        method of its base class to compute the spike-triggered average
        between the stimulus and response. Because raw voltage data is
        received, and not spike times, the class uses a "spike finder" to
        simulate spike times from the raw data.

        The method used for doing this depends on the `spike_method` parameter
        to the constructor, and either finds upward crossings of a threshold,
        or finds peak values above that threshold.
        '''
        super().__init__(stimfile, frame_rate, length, channel)

        # Determine spike-finder function
        self._spike_method = spike_method
        if spike_method == 'peaks':
            self._spike_finder = self._peak_finder
        elif spike_method == 'thresh':
            self._spike_finder = self._threshold_crossing_finder
        else:
            raise ValueError('''The spike-finding method {} is 
                    not supported.'''.format(spike_method))
        self.total_spike_count = 0


    def __str__(self):
        return 'OnlineSpikeTriggeredAverage(stimfile={}, length={:0.2f}, channels={})'.format(
                self._stimfile_name, self.length, self.channels)

        
    def update(self, dataframe):
        '''Update the online spike-triggered average.

        Parameters
        ----------

        dataframe : DataFrame object
            The chunk of data used to compute the STA.

        Raises
        ------

        If the time limits of the data frame are outside the bounds
        of the stimulus (e.g., if the recording ends after the stimulus),
        a ValueError is raised.

        Notes
        -----

        This method smooths the raw voltage data, uses the spike-time
        finder specified in the constructor to estimate "spike" times from
        the voltage data, and then updates the current result with the 
        spike-triggered average of the stimulus given those spike times.
        '''
        self.current_time = dataframe.stop

        # Smooth voltage data
        sample_rate = int(np.round(dataframe.nsamples() / 
                (dataframe.stop - dataframe.start)))
        data = dataframe[self.channels, :]
        smoothed = self._smooth_data(data, sample_rate)

        # Get frames for this chunk of data
        stim = self.frames_from_time(dataframe.start, dataframe.stop)
        stim_times = self.frame_times(dataframe.start, dataframe.stop)

        # Upsample stimulus to desired filter sample rate
        factor = int(np.round(self.filter_sample_rate / self.frame_rate))
        stim_upsampled, time_upsampled = pyret.stimulustools.upsample(
                stim, factor, stim_times)

        # Compute "spike" times
        spike_idx = self._spike_finder(smoothed)
        if len(spike_idx) == 0:
            return
        self.total_spike_count += len(spike_idx)
        tax = np.linspace(dataframe.start, dataframe.stop, dataframe.nsamples())
        spike_times = tax[spike_idx]

        # Compute STA
        new_sta = pyret.filtertools.sta(time_upsampled, stim_upsampled, 
                spike_times, self.filter_length)[0]
        if np.any(np.isnan(new_sta)):
            return
        self.current_result += new_sta[::-1, ...]
    
        # Decompose if necessary
        if self.stimulus.ndim == 3:
            self._spatial_rf, self._temporal_rf = pyret.filtertools.decompose(
                    self.current_result)

        self.current_title = '{0:#d}s ({1:#d} spikes)'.format(
                int(self.current_time), self.total_spike_count)


    def finalize(self):
        '''Finalize the online STA estimate by smoothing.'''
        self.current_result = np.squeeze(pyret.filtertools.smooth(
                np.atleast_3d(self.current_result)))
        if self.stimulus.ndim == 3:
            self._spatial_rf, self._temporal_rf = pyret.filtertools.decompose(
                    self.current_result)


    def _peak_finder(self, smoothed):
        '''Return indices of peaks in the single channel of data.'''
        threshold = self._compute_threshold(smoothed)
        diff = np.append(np.diff(smoothed), 0.0)
        peaks = np.append( (diff[:-1] < 0) & (diff[1:] > 0), 0)
        return np.where((smoothed >= threshold) & peaks)[0]


    def _threshold_crossing_finder(self, smoothed):
        '''Return indices of upward threshold crossings in the single
        channel of data.
        '''
        threshold = self._compute_threshold(smoothed)
        above = smoothed > threshold
        diff = np.append(np.diff(above), 0)
        upward_crossings = np.where((diff[:-1] < 0) & (diff[1:] > 0))[0]
        return upward_crossings


    def _compute_threshold(self, data, percentile=95):
        '''Return a simple threshold for spike-detection.'''
        return np.median(np.abs(data)) * 4.5


    def _smooth_data(self, data, sample_rate, lp_cutoff=10, hp_cutoff=None):
        '''Smooth raw voltage data, by band-pass filtering.

        Parameters
        ----------

        data : array_like
            The data to be smoothed.
    
        sample_rate : float
            The sample rate of the data.

        lp_cutoff, hp_cutoff : floats
            The low and high-pass cutoff frequencies for the Butterworth filter.

        Returns
        -------

        smoothed : array_like
            The band-pass filtered data.
        '''

        fs = sample_rate // 2
        if hp_cutoff is None:
            hp_cutoff = fs - 0.1 * fs
        order = 6
        frequencies = np.array((lp_cutoff, hp_cutoff)) / fs
        b, a = signal.butter(order, frequencies, btype='band')
        return signal.filtfilt(b, a, data)

