from __future__ import division  # so that 1/3=0.333 instead of 1/3=0
from psychopy import locale_setup, visual, core, data, event, logging, sound, gui
from psychopy.constants import *  # things like STARTED, FINISHED
import numpy as np  # whole numpy lib is available, prepend 'np.'
from numpy import sin, cos, tan, log, log10, pi, average, sqrt, std, deg2rad, rad2deg, linspace, asarray
from numpy.random import random, randint, normal, shuffle
import math # math needed for sqrt
from reach_constants_Su16 import *
#some useful functions
def getPokeType(penDist, targDist, tooSoonBool):
	#determines whether poke occured within target region, penalty region, overlap region, or outside of stimulus structure
    if (penDist > PEN_TARG_RADIUS) and (targDist < PEN_TARG_RADIUS):
    	poke_type = TARG_POKE
    elif (penDist < PEN_TARG_RADIUS) and (targDist > PEN_TARG_RADIUS):
    	poke_type = PEN_POKE
    elif (penDist < PEN_TARG_RADIUS) and (targDist < PEN_TARG_RADIUS):
    	poke_type = MIX_POKE
    elif (penDist > PEN_TARG_RADIUS) and (targDist > PEN_TARG_RADIUS):
    	poke_type = MISS_POKE
    elif tooSoonBool:
    	poke_type = EARLY_POKE
    else:
    	poke_type = LATE_POKE
    return poke_type 

def getTrialScore(pokeType, rewardVal, penaltyVal):
    #assign points based on what region the poke occured within
    if pokeType == EARLY_POKE:
        score = float('nan')
    elif pokeType == TARG_POKE:
        score = rewardVal # MAGIC NUMBERS  
    elif pokeType == PEN_POKE:
        score = penaltyVal
    elif pokeType == MIX_POKE:
        score = rewardVal + penaltyVal
    elif pokeType == MISS_POKE:
        score = 0
    else: 
        score = TOO_SLOW_PEN
    return score

def getTrialAcc(pokeType):
    #assign accuracy based on what region the poke occured within
    if pokeType == TARG_POKE:
        accuracy = 1
    elif pokeType == PEN_POKE:
        accuracy = -1
    elif pokeType == MIX_POKE:
        accuracy = 0
    elif pokeType == MISS_POKE:
        accuracy = -2
    else: 
        accuracy = float('nan')
    return accuracy

def getOpacities(pokeType, frameCount):
    if pokeType == MIX_POKE and not (frameCount % 2):
        targ_opac = 1
        pen_opac  = 1
    elif pokeType == MIX_POKE  and (frameCount % 2): 
        targ_opac = 0
        pen_opac  = 0
    elif pokeType == PEN_POKE  and not (frameCount % 2):
        targ_opac = 1
        pen_opac  = 1
    elif pokeType == PEN_POKE and (frameCount % 2):
        targ_opac = 1
        pen_opac  = 0
    elif pokeType == TARG_POKE and not (frameCount % 2):
        targ_opac = 1
        pen_opac  = 1
    elif pokeType == TARG_POKE and (frameCount % 2):
        targ_opac = 0
        pen_opac  = 1
    elif pokeType == MISS_POKE:
        targ_opac = 1
        pen_opac  = 1
    else:
    	targ_opac = 1
    	pen_opac  = 1
    return [targ_opac, pen_opac]
    
def dist(p1, p2):
	return math.sqrt((p2[0] - p1[0]) ** 2 + (p2[1] - p1[1]) ** 2)