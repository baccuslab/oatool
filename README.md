# `oatool`

Library and simple tool for online analysis of retinal data in the Baccus Lab.

(C) 2017 Benjamin Naecker bnaecker@stanford.edu

## Overview

Online analysis is the process of analyzing data from retinal experiments
in roughly real-time. The most common example is an online receptive field
mapping of a cell recorded intracellularly. This may tell the experimenter
what type of cell it is, and thus to continue or abandon the recording, or
to show some specific stimulus, etc.

`oalib` consists of a small Python library, which implements a few of the
most basic analyses, and an executable for using those analyses from the
command-line and plotting the results.

## Dependencies

`oatool` is designed to work with the Baccus lab recording software pipeline.
This is installed already on the Fairchild D239 recording computers, but you
can easily install and build it on your own machine for testing. The quickest
way to do this is to install the recording software super-repository from
the lab GitHub site: 
[recording-software](https://github.com/baccuslab/recording-software)

The pipeline has dependencies of its own, see that repository for information
about how to get those. On Linux, all of them should be installable via
the system's package manager (`apt` or similar), and you can use 
[Homebrew](https://brew.sh) to install everything on macOS.

The Python dependencies are as follows:

- [Python3](https://python.org)
- [`pyret`](https://pyret.readthedocs.io) (which has several of its own)
- [`matplotlib`](https://www.matplotlib.org)
- [`h5py`](http://www.h5py.org)
- [`bldsclient`](https://github.com/baccuslab/libblds-client) 
	(part of the `libblds-client` Baccus lab private repository)

### Python path

The required repositories must be in your Python module search path. Most of
the packages can be installed via `pip`, which will put them on the path
by default. The exception is the `bldsclient` module, which is only available
by downloading from GitHub. You can add it to your path in three ways:

- `cd /path/to/libblds-client/bldsclient/` and `pip3 install .`
- importing the `sys` module in your scripts, and then doing:
	`sys.path.append('/path/to/libblds-client/python')`
- updating the shell environment variable `PYTHONPATH` by doing:
	`PYTHONPATH="$PYTHONPATH:/path/to/libblds-client/python"`

## `oatool`

`oatool` is what most people will use. It is a simple executable with which
experimenters can run the currently-implemented online analyses in live
experiments.

### Usage

	oatool.py [ -h | --help]
		--stimfile FILE
		[--frame-rate RATE] [--hostname HOST]
		[--analysis-type TYPE] [--spike-finder TYPE]
		[--duration DURATION] [--analysis-length LEN]
		[--interval INTERVAL] [--custom-analysis ANALYSIS_MODULE]
		[--no-wait]

The most salient argument is `--stimfile`, which gives the stimulus
currently being used in the experiment. This must be an HDF5 file, and
contain a dataset called `'stimulus'`. If the frame rate is not specified
on the command line via the `--frame-rate` parameter, the `'stimulus'`
dataset must also have an attribute called `'frame-rate'` which does specify
the rate in Hz.

The other important argument is the `--analysis-type`, which selects
which online analysis to run. The two current choices are `'revcorr'` for
an online reverse-correlation of a membrane potential with the stimulus, and
`'sta'` for an online spike-triggered average. Others may be implemented
in the future.

## `oalib`

`oalib` contains an actual Python module which implements all the online
analysis code. The module is designed to be flexible and easily extensible.

Analyses are actually Python classes, which must inherit from the abstract
base class `OnlineAnalysis`. This defines a basic structure of an analysis,
and has two abstract methods called, `update()` and `plot()`. The `update`
method must take a frame of data and use that to update its running analysis.
What that actually means is totally up to the subclass. The `plot` method
should actually plot the results.

There is another abstract class `OnlineReceptiveField`, which specializes
the `OnlineAnalysis` base class for implementing online receptive field
computation. This really just implements the `plot` method, to appropriately
plot 1, 2, or 3D receptive fields as they are computed.

The `OnlineReverseCorrelation` and `OnlineSpikeTriggeredAverage` classes
implement the `'revcorr'` and `'sta'` analyses that can be currently selected
by the `--analysis-type` flag to `oatool`.

## Extending

It should be pretty straightforward to write a new analysis. Create a new class
inherited from `OnlineAnalysis` (or `OnlineReceptiveField` if you want), and
implement the `update` and `plot` methods. There are more details in the
docstrings for those classes. But here's a quick example of how to do this,
to create an analysis that just estimates the total activity of every channel.

	class OnlineActivityMap(OnlineAnalysis):
		def __init__(self, stimfile, frame_rate, length, channels):
			super().__init__(stimfile, frame_rate, length, channels)
			self.current_result = np.zeros((len(channels),))

			# We'll plot the activity over time as a bar graph
			self.axes = self.figure.add_subplot(111)
			self.activity_map = self.axes.bar(self.channels, self.current_result)

		def update(self, dataframe):
			# "activity" is just the cumulative variance of the raw data
			self.current_result += dataframe.data.var(axis=1)

		def plot(self):
			for idx, bar in enumerate(self.activity_map):
				bar.set_height(self.current_result[idx])
			plt.pause(0.01) # Force update of graphics stack

This hasn't been tested, and most likely won't work. But it does give all
the necessary details for how to actually include a new online analysis.

The analysis you write can either be included directly in `oalib`, or can
be in your own module. You can then give the path to the module as the
argument to the `--custom-analysis` parameter to `oatool`. This will verify
that there is a class in that file that inherits from `OnlineAnalysis`,
and will use that instead. This selection mechanism is pretty dumb, it just
chooses the first class that has the right parentage. So only write one
per file.

