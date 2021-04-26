# This file will contain all the constants and components for the Reach Experiment
from __future__ import division  # so that 1/3=0.333 instead of 1/3=0
from psychopy import locale_setup, visual, core, data, event, logging, sound, gui
from psychopy.constants import *  # things like STARTED, FINISHED
import numpy as np  # whole numpy lib is available, prepend 'np.'
from numpy import sin, cos, tan, log, log10, pi, average, sqrt, std, deg2rad, rad2deg, linspace, asarray
from numpy.random import random, randint, normal, shuffle
import os  # handy system and path functions
import sys # to get file system encoding
import pygame
import math # math needed for sqrt
import pandas as pd
from pandas import DataFrame, Series
import random
from pyglet.window import key, Window
import winsound


# Store info about the experiment session
expName = u'malohney'  # from the Builder filename that created this script
expInfo = {u'session': u'001', u'participant': u''}
dlg     = gui.DlgFromDict(dictionary=expInfo, title=expName)
if dlg.OK == False: core.quit()  # user pressed cancel
expInfo['date']    = data.getDateStr()  # add a simple timestamp
expInfo['expName'] = expName

# Ensure that relative paths start from the same directory as this script
_thisDir = os.path.dirname(os.path.abspath(__file__)).decode(sys.getfilesystemencoding())
os.chdir(_thisDir)


#this reads in the data frame of circle configurations
testInfo   = pd.DataFrame.from_csv('testTrials.csv')
train2Info = pd.DataFrame.from_csv('trainTrialsPt2.csv')

#set standard timeLimit. This assumes it is the same for every trial in block
TRIAL_TIME_LIMIT = testInfo.iloc[0]['timeLimit']

#Use this to turn on and off the tone signal
USE_BEEP = False

# Data file name stem = absolute path + name; later add .psyexp, .csv, .log, etc
filename = _thisDir + os.sep + u'data/%s_%s_%s' %(expInfo['participant'], expName, expInfo['date'])

endExpNow = False  # flag for 'escape' or other condition => quit the exp

#durations
TRAIN_INSTR_DUR   = 5
TRAIN_FEED_DUR    = 1
TEST_INSTR_DUR    = 7
TRAIN2_INSTR_DUR  = 7
VALUE_DISPLAY_DUR = 1.5
TRAIN2_FEED_DUR   = 1
TEST_FEED_DUR     = 1
POINT_FEED_DUR    = 3
#PAYMENT_FEED_DUR  = 10


#loop sizes
TOTAL_TRAIN_BLOCK       = 8 #8
TRAIN_TRIALS_PER_BLOCK  = 36 #36

TOTAL_TRAIN2_BLOCK      = 8 #8
TRAIN2_TRIALS_PER_BLOCK = 36 #36


TOTAL_TEST_BLOCK       = 10 #This number can be less than the number of blocks in dataframe but not more
TEST_TRIALS_PER_BLOCK  = 60 #int(len(testInfo.index) / (testInfo.iloc[len(testInfo.index) - 1]["blockNum"] + 1)) #length of dataframe / number of blocks in dataframe

#RGB Colors
WHITE  = [1, 1, 1]
BLUE   = [-1, -1, 1]
GREEN  = [-1, 1, -1]
RED    = [1, -1, -1]
YELLOW = [1, 1, -1]
PURPLE = [1, -1, 1]
GREY   = [0, 0, 0]
GOLD   = [1, .84, 0]
ORANGE = [1, .5, 0]
LBLUE  = [0, 1, 1]

#radius of the circles
PEN_TARG_RADIUS = 32

#timeout penalty
TOO_SLOW_PEN = -5

#coordinates for bottom center of window
BOT_CENTER   = [0, -899]

#poke_type codes
TARG_POKE  = 1
PEN_POKE   = 2
MIX_POKE   = 3
MISS_POKE  = 4
LATE_POKE  = 5
EARLY_POKE = 0

#create lists for storing training data
aimDotPositions     = []
pokePositionsTrain  = []
pokePositionstrain2 = []
pokeDistancesTrain  = []
trialNumsTrain      = []
pokeTimesTrain      = []
pokeTimestrain2     = []            
releaseTimesTrain   = []
releaseTimesTrain2  = []

#create lists to store test data
scores            = []
too_slow_bools    = []
totalScores       = []
pokePositionsTest = []
penaltyPositions  = []
targetPositions   = []
penaltyDistances  = []
targetDistances   = []
pokeTimesTest     = [] #rename to have RT
trialNumber       = []
blockNumber       = []
targetSizes       = []
sepDistances      = []
rewardValue       = []
penaltyValue      = []
timeLimits        = []
withinStimBools   = []  # 1 == poked within stimulus config. , 0 == poked outside stim. config. , NaN if no poke
accuracies        = []  # 1 == target poke, 0 == middle poke, -1 == penalty poke, NaN == else (i.e., too slow or outside config.)
tooSlowBools      = []  # 1 == too slow, 0 == not too slow
releaseTimesTest  = []
pokeTypes         = []