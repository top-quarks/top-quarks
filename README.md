Hello!

Below you can find a outline of how to reproduce my solution for the TrackML competition.
If you run into any trouble with the setup/code or have any questions please contact me at top-quarks@protonmail.com

#HARDWARE: (The following specs were used for testing and creating most submissions*)
Ubuntu 16.04 LTS
Lenovo Legion Y720 (4 core i7-7770HQ, 16 GB memory, didn't use the GTX 1060 graphics card)

* For all submissions and test except a few of the last ones I used this hardware. For creating the submission I uploaded I actually used a remote computer, but there was no reason to do so, and I re-ran it on this hardware to make sure there was no difference. So I'll just present numbers from this hardware.

The minimum hardware is currently limited by the 4.07GB peak memory usage of the largest event. The total run-time for all 125 events is assumed to be approximately 16 hours divided by the number of cores.

#SOFTWARE
The final model is standalone C++ requiring only standard STL libraries and C++11 support, python is only used for scripting (re-training, running multiple events in parallel, stitching submission files)
g++ (Ubuntu 5.5.0-12ubuntu1~16.04) 5.5.0 20171010
Python 3.5.2

For optional re-training of logistic regression parameters
scikit-learn==0.19.1
numpy==1.14.3

#DATA SETUP
The model assumes it can find the following (un-zipped) competition files (the contents of blacklist_training.zip is contained in blacklist/)

blacklist/
detectors.csv
train_100_events/
test/

in a folder specified by "base_path" in input.hpp . There is no fatal error if blacklist isn't found, and we don't need train_100_events/ except for re-training the model. Also, note that "base_path" can't use ~/ to point to home, it needs a full or relative path like "/home/icecuber/kaggle_data" or "../kaggle_data/".

#MODEL RUN
The model is run with

python3 predict.py

This took all little under 4 hours on my hardware, running 4 events in parallel at any time (you may adjust the "parallel" setting in run_test.py. It never came close to the 16 GB memory available. The largest event had a peak usage of 4.07 GB, the average event had peak usage 2.78 GB.

The output (which is overwritten every time) is a file called "submissions/submission.csv" which could be uploaded to Kaggle (if gzip1.6 is available, there will also be a compressed version "submissions/submission.csv.gz" which was the one I uploaded. Otherwise you might get an unimportant error message saying it couldn't compress it).

The first lines of run_test.py are (hopefully) very intuitive, so take a look if you want to f.ex. run on a subset of the events.

The logistic regression model for pairs is re-trained with

python3 train.py

Expect a lot of ugly output, rm fails are expected. train.py will overwrite the old training in trained/pair_logreg3