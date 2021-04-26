3# Loss aversion task
# Displays a series of gambles that participants can either accept or reject

import random, pygame, sys, time
import pandas as pd
import numpy as np
from pandas import Series, DataFrame
from pygame.locals import *

#             R    G    B
WHITE     = (255, 255, 255)
BLACK     = (  0,   0,   0)
RED       = (255,   0,   0)
GREEN     = (  0, 255,   0)
DARKGREEN = (  0, 155,   0)
DARKGRAY  = ( 40,  40,  40)
BGCOLOR = BLACK

WINDOWWIDTH = 1500
WINDOWHEIGHT = 900
CELLSIZE = 2
RESPONSEKEYFONT = 50
OFFERFONT = 175
assert WINDOWWIDTH % CELLSIZE == 0, "Window width must be a multiple of cell size."
assert WINDOWHEIGHT % CELLSIZE == 0, "Window height must be a multiple of cell size."
CELLWIDTH = int(WINDOWWIDTH / CELLSIZE)
CELLHEIGHT = int(WINDOWHEIGHT / CELLSIZE)

LOSSES = range(5,16) 
GAINS = [x*2 for x in LOSSES]
    
def main():
    global BASICFONT, DISPLAYSURF
    
    subject = raw_input("Please enter your subject number:")
    #initialize pygame
    pygame.init()
    DISPLAYSURF = pygame.display.set_mode((WINDOWWIDTH, WINDOWHEIGHT))
    BASICFONT = pygame.font.Font('freesansbold.ttf', 18)
    
    
    #shuffle losses and gains into random orders
    shuffledLosses = shuffle_array(LOSSES)
    shuffledGains = shuffle_array(GAINS)
    
    allOffers = []
    responses = []
    acceptedGambles = []
    rejectedGambles = []
    offerResponses = []
    
    # create an array of paired gains and losses in every combination
    for gind,gain in enumerate(shuffledGains):
        for lind,loss in enumerate(shuffledLosses):
            allOffers.append([gain, loss])
    shuffledOffers = shuffle_array(allOffers) # shuffle the order of the offers
    
    for oind, offer in enumerate(shuffledOffers):
        
        DISPLAYSURF.fill(BGCOLOR)
        showResponseKeys()
        pygame.display.update()
        pygame.time.wait(500)
        showMessage(makeTextOffer(offer[0], offer[1]))
        response ='none'
        while response == 'none':

            for event in pygame.event.get(): # event handling loop
                if event.type == QUIT:
                    terminate()
                elif event.type == KEYDOWN:
                    response = checkKeyPress(event.key)
        responses.append(response)
        if response == 'accept':
            acceptedGambles.append(offer)
            offerResponses.append(1)
        else:
            rejectedGambles.append(offer)
            offerResponses.append(0)
    pygame.quit()
    
    # pick one trial at random that the participant accepted
    randtrial = random.randint(0, len(offerResponses)-1)
    if offerResponses[randtrial] == 1:
        print 'The randomly selected offer you accepted was $%s' % \
            makeTextOffer(shuffledOffers[randtrial][0],shuffledOffers[randtrial][1])
        HorT = raw_input("Please call heads or tails\n")
        while HorT != 'heads' and HorT != 'tails':
            HorT = raw_input("Please type 'heads' or 'tails'\n")
        choice = random.choice(['heads', 'tails'])
        print 'flipping coin . . .'
        time.sleep(2)
        
        print choice
        if HorT == choice:
            print 'Congratulations! You won $%i' % shuffledOffers[randtrial][0]
        else:
            print 'Unfortunately you lost $%i' % shuffledOffers[randtrial][1]
    else:
        print 'You rejected the randomly selected offer: %s' % \
            makeTextOffer(shuffledOffers[randtrial][0],shuffledOffers[randtrial][1])
    
  
    np.save(subject + "_accepted.npy", acceptedGambles)
    np.save(subject + "_rejected.npy", rejectedGambles)
    terminate()
        
        
    
def checkKeyPress(key):
    #keytype = 'error'
    if (key == K_1 or key == K_y):
        keytype = 'accept'
    elif (key == K_0 or key == K_n):
        keytype = 'reject'
    elif key == K_ESCAPE:
        terminate()
    else:
        keytype = 'none'
    return keytype

def checkForKeyPress():
    if len(pygame.event.get(QUIT)) > 0:
        terminate()

    keyUpEvents = pygame.event.get(KEYUP)
    if len(keyUpEvents) == 0:
        return None
    if keyUpEvents[0].key == K_ESCAPE:
        terminate()
    return keyUpEvents[0].key
    
    
def terminate():
    pygame.quit()
    sys.exit()

def drawTrialNumber(trialNum):
    trialSurf = BASICFONT.render('Trial # %s' % (trialNum), True, WHITE)
    trialRect = trialSurf.get_rect()
    trialRect.topleft = (20, 10)
    DISPLAYSURF.blit(trialSurf, trialRect)

def makeTextOffer(gain, loss):
    textOffer = '+$' + str(gain) +' | ' + '-$' + str(loss)
    return textOffer
    
def showMessage(message):
    gameOverFont = pygame.font.Font('freesansbold.ttf', OFFERFONT)
    gameSurf = gameOverFont.render(message, True, WHITE)  
    gameRect = gameSurf.get_rect()
    gameRect.midtop = (WINDOWWIDTH / 2, WINDOWHEIGHT / 3)

    DISPLAYSURF.blit(gameSurf, gameRect)
    #drawPressKeyMsg()
    pygame.display.update()
    pygame.time.wait(500)
    checkForKeyPress() # clear out any key presses in the event queue
    
    #~ while True:
        #~ if checkForKeyPress():
            #~ pygame.event.get() # clear event queue
            #~ return
def showResponseKeys():
    keyFont = pygame.font.Font('freesansbold.ttf', RESPONSEKEYFONT)
    keySurf = keyFont.render('Press 1 to accept or 0 to reject', True, WHITE)
    keyRect = keySurf.get_rect()
    keyRect.midtop = (WINDOWWIDTH/2, WINDOWHEIGHT*5/6)
    DISPLAYSURF.blit(keySurf, keyRect)
    pygame.display.update()
    
def shuffle_array(array):
    random.shuffle(array)
    return array

if __name__ == '__main__':
    main()