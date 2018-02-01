#!/usr/bin/env python3
'''oatool.py
Tool for online analysis of retinal data.
(C) 2017 Benjamin Naecker bnaecker@stanford.edu
'''

import argparse
import sys
import time
import warnings

import h5py

import bldsclient
import oalib

def parse_command_line():
    '''Parse command-line arguments to the oatool.'''
    parser = argparse.ArgumentParser(description='''
            This tool collects data in chunks, and using the provided stimulus
            file, runs one of several online analyses. The analysis will be
            plotted as new data is collected.

            Users must specify a stimulus file, but all other arguments are
            optional. If no analysis type is specified, the tool will assume
            an online intracellular receptive field is desired. This will be
            computed via reverse-correlation between the stimulus and response.''')

    parser.add_argument('--stimfile', 
            help='''HDF5 file containing the stimulus against 
            which the analysis will be run. This argument is required.''',
            required=True)
    parser.add_argument('--frame-rate',
            help='''Effective frame rate of the stimulus, in Hz.
            If a value is not given, it must be available as an
            attribute of the "stimulus" dataset in the given 
            stimulus file.''',
            default=0.0,
            type=float)
    parser.add_argument('--hostname',
            help='''Hostname or IP address of machine running BLDS.
            Defaults to "localhost".''',
            default='localhost')
    parser.add_argument('--analysis-type',
            help='''Type of analysis to run. Built-in options are:
            "revcorr", for an intracellular receptive field via 
            reverse-correlation, and "sta", for an extracellular
            receptive field via spike-triggered averaging.''',
            choices=('revcorr', 'sta'),
            default='revcorr')
    parser.add_argument('--spike-finder',
            help='''Method for finding "spikes" in the data, when using the "sta"
            analysis method. There are two implemented methods. A value of "thresh"
            will use a simple simple upward threshold crossing to determine 
            spike times. This is fast, but may be less precise. A value of "peaks"
            will smooth the data and find the times of peaks. This is slower, but
            will be more accurate. Note that this parameter has no effect when the
            "revcorr" analysis type is used.''',
            choices=('thresh', 'peaks'),
            default='peaks')
    parser.add_argument('--channels',
            help='''Channels in the dataset on which the analysis will
            be run. This should be specified as a comma-, or dash-separated 
            list of channels, e.g., 0,1,2 or 0-8,12-14. If not given, defaults 
            to channel 1, which is the intracellular voltage channel in the MCS
            recording system. The special value of "all" can be passed to indicate
            that all data channels (however there happen to be) should be passed
            to the analysis routines.''',
            default='1')
    parser.add_argument('--duration',
            help='''Total amount of data to analyze, in seconds. If not given,
            data will be collected until either the analysis is interrupted or
            the recording ends.''',
            default=-1.0,
            type=float)
    parser.add_argument('--analysis-length',
            help='''Length of the analysis, in seconds. The exact meaning of
            this parameter depends on the analysis type. For most, it is the 
            length of an online filter to estimate. But for custom analyses,
            this may be used in other ways.
            ''',
            default=0.50,
            type=float)
    parser.add_argument('--interval',
            help='''Interval at which data is collected and the analysis is
            updated. This must be at least twice as long as the length of the
            analysis, which is the default. Any value less than this is set
            to that length.
            ''',
            default=0.0,
            type=float)
    parser.add_argument('--custom-analysis',
            help='''Name of a Python file specifying a custom online analysis
            to be run. This module must contain an online analysis subclass,
            which implements the appropriate methods. See the "oalib" module 
            for more details.''')
    parser.add_argument('--no-wait',
            help='''By default, the tool will wait until a recording is started
            before collecting data and analyzing it. Passing this flag instructs
            the tool to actually start the recording.''',
            action='store_true', 
            default=False)

    ns = parser.parse_args(sys.argv[1:])
    chanlist = parse_channel_list(ns.channels)

    return (ns.stimfile, ns.frame_rate, ns.hostname,
            ns.analysis_type, chanlist, ns.duration, ns.analysis_length, 
            ns.interval, ns.spike_finder, ns.custom_analysis, ns.no_wait)


def parse_channel_list(channels):
    '''Parse a comma- and/or dash-separated list of channels into 
    a tuple of integers.
    '''
    chanlist = set()
    groups = channels.split(',')
    for group in groups:
        chans = group.split('-')
        if len(chans) == 1:
            chanlist.add(int(chans[0]))
        elif len(chans) == 2:
            chanlist.update(set(range(*map(int, chans))))
        if len(chans) > 2:
            raise ValueError('''The given channel list ({}) 
                    could not be parsed. Specify the list as a comma- and dash-
                    separated list of channels.'''.format(channels))
    return tuple(chanlist)


def connect_and_verify_source(hostname, chanlist, length):
    '''Connect to the Baccus lab data server, and verify that the requested
    channel and analysis duration are valid. Note that this will block until
    a data source has been created and is ready to record.
    '''
    # Connect client
    client = bldsclient.BldsClient(hostname)
    client.connect()

    # Wait for data source to become available
    while not client.get('source-exists'):
        time.sleep(1.0)
   
    # Verify the channels are valid
    nchannels = client.get_source('nchannels')
    for channel in chanlist:
        if (channel < 0) or (channel >= nchannels):
            raise ValueError('''The requested analysis channel is not 
                    valid for the current source.''')

    # Verify the requested length is valid, or truncate to valid
    max_length = client.get('recording-length')
    if length < 0.0:
        length = max_length
    elif length > max_length:
        warnings.warn('''Requested analysis length is longer than recording,
                truncating to {:0.2f} seconds'''.format(max_length))
        length = max_length
    return client, length


def create_online_analysis(analysis_type, custom_module,
        stimfile, frame_rate, chanlist, analysis_length, spike_finder):
    '''Create an online analysis object of the requsted type.'''
    oa = None
    if custom_module is None:
        if analysis_type == 'revcorr':
            oa = oalib.OnlineReverseCorrelation(stimfile, frame_rate, 
                    analysis_length, chanlist)
        elif analysis_type == 'sta':
            oa = oalib.OnlineSpikeTriggeredAverage(stimfile, frame_rate,
                    analysis_length, chanlist, spike_finder)
        else:
            raise NotImplementedError('''The requested analysis type '{}'
                    is not implemented.'''.format(analysis_type))
    else:
        # Find the OnlineAnalysis subclass implemented in the requested
        CustomClass = find_custom_analysis_class(custom_module)
        oa = CustomClass(stimfile, frame_rate, analysis_length, chanlist)

    return oa


def find_custom_analysis_class(pyfile):
    '''Try to find a custom analysis class in the given module.

    This allows users to specify their own Python module with a new
    analysis class, which can then use the existing oatool infrastructure.

    Parameters
    ----------

    pyfile : str
        The absolute or relative path of a Python module containing
        the analysis class to be run.

    Raises
    ------

    A FileNotFoundError is raised if the module couldn't be found.
    A TypeError is raised if the module exists but doesn't appear
    to contain a class or that class does not inherit from the abstract
    OnlineAnalysis base class.
    '''
    import inspect
    import importlib.util
    import os.path

    if not os.path.exists(pyfile):
        raise FileNotFoundError('The given module does not exist')

    modname = inspect.getmodulename(pyfile)
    if modname is None:
        raise ValueError('The Python file is not a module')
    spec = importlib.util.spec_from_file_location(modname, pyfile)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    # Search for classes which inherit from OnlineAnalysis
    for name, cls in inspect.getmembers(module, inspect.isclass):
        if issubclass(cls, oalib.OnlineAnalysis):
            return cls

    raise ValueError('''The module {} does not contain a subclass of
            oalib.OnlineAnalysis'''.format(module.__name__))

if __name__ == '__main__':

    warnings.simplefilter('ignore')

    # Parse command-line arguments
    (stimfile, frame_rate, hostname, analysis_type, 
            chanlist, duration, analysis_length, interval, spike_finder,
            custom_analysis, no_wait) = parse_command_line()

    try:

        # Connect to the data server and verify the requested channels and duration
        client, duration = connect_and_verify_source(hostname, chanlist, duration)

        # Create analysis object of requested type
        oa = create_online_analysis(analysis_type, custom_analysis,
                stimfile, frame_rate, chanlist, analysis_length, spike_finder)
        print('Running online analysis for {:#0.2f} seconds using: {}'.format(
                duration, oa))

        # Wait for the start of the recording
        if no_wait:
            if not client.get('recording-exists'):
                try:
                    client.start_recording()
                except bldsclient.BldsError:
                    pass
                print('Starting recording')
        else:
            print('Waiting for recording to start')
            while not client.get('recording-exists'):
                time.sleep(1.0)

        # Verify duration and interval
        duration = min(duration, client.get('recording-length'))
        interval = max(interval, 2 * analysis_length)

        # Collect data and run the online analysis
        print('Starting data stream and online analysis')
        position = 0.0
        while position < duration:
            frame = client.get_data(position, position + interval)
            oa.update(frame)
            oa.plot()
            position += interval

        # Notify finished without error
        print('Analysis and/or recording complete.')

    except KeyboardInterrupt:
        print('Analysis has been stopped.')

    except ValueError as exc:
        print('Analysis ended early, most likely because the recording stopped.')
        print(exc)

    except Exception as exc:
        import traceback
        print('An unexpected error occurred: {}'.format(exc))
        print(traceback.print_exc())

    finally:

        # Finalize the analysis
        oa.finalize()
        oa.plot()

        # Cleanup
        client.disconnect()

        print('Press ENTER to quit')
        sys.stdin.read(1)

