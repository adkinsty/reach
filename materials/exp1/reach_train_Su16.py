#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
This experiment was created using PsychoPy2 Experiment Builder (v1.83.04), May 11, 2016, at 14:16
If you publish work using this script please cite the relevant PsychoPy publications
  Peirce, JW (2007) PsychoPy - Psychophysics software in Python. Journal of Neuroscience Methods, 162(1-2), 8-13.
  Peirce, JW (2009) Generating stimuli for neuroscience using PsychoPy. Frontiers in Neuroinformatics, 2:10. doi: 10.3389/neuro.11.010.2008
"""

from reach_constants_Su16 import *
from reach_functions_Su16 import *
from reach_components_Su16 import *

# Ensure that relative paths start from the same directory as this script
_thisDir = os.path.dirname(os.path.abspath(__file__)).decode(sys.getfilesystemencoding())
os.chdir(_thisDir)

current_train_block = 0

mouse.setVisible(0)

#Loop for part 1 of training (single target); each iteration is one block
while current_train_block < TOTAL_TRAIN_BLOCK:

    current_train_trial = int(TRAIN_TRIALS_PER_BLOCK * current_train_block)
    total_train_trial   = current_train_trial + TRAIN_TRIALS_PER_BLOCK
    poke_num = 0

    train_instr_timer.reset()
    train_instr_timer.add(TRAIN_INSTR_DUR)
    train_instr_text = 'You are about to begin block ' + str(current_train_block + 1) + ' of ' + str(TOTAL_TRAIN_BLOCK) + '.\nWhen the yellow dot appears, try to poke it with your finger. Please press Z to continue.'
    train_instructions.setText(train_instr_text)

    while not keyboard.state.has_key('z'):
        
        train_instructions.draw() 
        win.flip()

        # check for quit (the Esc key)
        if keyboard.state.has_key('escape'):
            win.close()
            core.quit()

    #Loop for each block in part 1 of training
    while current_train_trial < total_train_trial:

        poke_now          = False
        poked             = False
        released_too_soon = False
        
        #------Prepare to start Routine "training"-------
        t = 0
        trainingClock.reset()  # clock 
        frameN = -1

        #random or uniform distribution?
        targ_pos_x = testInfo.iloc[current_train_trial]['x_position']
        targ_pos_y = testInfo.iloc[current_train_trial]['y_position']

        target.pos        = [targ_pos_x, targ_pos_y]
        pen.pos           = [targ_pos_x - PEN_TARG_RADIUS, targ_pos_y]
        train_aim_dot.pos = target.pos
        

        target.linecolor  = GREEN
        continueRoutine   = True
        space_pressed     = False

        #change fixation color to WHITE
        cross_vertical.fillColor   = WHITE
        cross_horizontal.fillColor = WHITE
        cross_vertical.lineColor   = WHITE
        cross_horizontal.lineColor = WHITE
        space_bar_timer.add(6)
        #FIXATION CROSS
        #Ends when spacebar is pressed
        while not keyboard.state.has_key(' '):
            # train2_text.setText(str(keyboard.state))
            # train2_text.draw()
            cross_vertical.draw()
            cross_horizontal.draw()
            press_text.draw()
            
            if space_bar_timer.getTime > 0:
                press_text.draw()
            else:
                please_press.draw()
                
            win.flip()

            # check for quit (the Esc key)
            if keyboard.state.has_key('escape'):
                win.close()
                core.quit()

        press_time = trainingClock.getTime()
        blue_fix_timer.reset()
        blue_fix_timer.add(random.uniform(.4, .6))

        #make fixation cross blue
        cross_vertical.fillColor = BLUE
        cross_horizontal.fillColor = BLUE
        cross_vertical.lineColor = BLUE
        cross_horizontal.lineColor = BLUE


        #BLUE FIXATION and BOUND BOX; Duration between 400 and 600 ms
        while blue_fix_timer.getTime() > 0:

            cross_vertical.draw()
            cross_horizontal.draw()
            bound_box.draw()
            win.flip()
            
            #the space bar should be depressed during this time; if not, the release is considered too soon
            if not keyboard.state.has_key(' '):
                released_too_soon = True

            # check for quit (the Esc key)
            if keyboard.state.has_key('escape'):
                win.close()
                core.quit()
              

        #initialize mouse position to bottom center of screen
        mouse.setPos(BOT_CENTER)
        aim_dot_timer.reset()
        aim_dot_timer.add(TRIAL_TIME_LIMIT)
        
        #bool for whether the time of release has been stored
        release_marked = False
        release_time_train = 0
        too_slow = False

        #Target presentation
        while not released_too_soon and not poked:

            # get current time
            t = TRIAL_TIME_LIMIT - aim_dot_timer.getTime() 
            frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
            if aim_dot_timer.getTime() <= 0:
                too_slow = True
            # check for quit (the Esc key)
            if keyboard.state.has_key('escape'):
                win.close()
                core.quit()
            # poke_now =
            # # true iff poke occured at this frame within boundaries of 'background', i.e., anywhere in window
            
            
            released = not keyboard.state.has_key(' ')
            if released and not release_marked:
                release_time_train   = t
                release_marked       = True
            #too soon condition, 100ms
            if release_time_train < .1 and release_time_train != 0:
                released_too_soon = True
                break

            #if poke occured at this frmae, then prepare to draw poke_dot at that the poke location at next frame and store time elapsed since trial start
            if  mouse.isPressedIn(background, buttons=[0]) and  not keyboard.state.has_key(' '): 
                poked            = True
                poke_pos         = mouse.getPos()
                poke_dot.pos     = poke_pos
                poke_time_train  = t # rename to rt of some sort
                break

            target.draw()    
            pen.draw()   
            train_aim_dot.draw()
            background.draw()
            win.flip()
            

        feedback_timer.reset()
        feedback_timer.add(TRAIN_FEED_DUR)
        
        #Post-poke feedback presentation
        while feedback_timer.getTime() > 0:

            if too_slow:
                no_poke_text.draw()  
                win.flip() 
            elif released_too_soon:
                too_soon_text.draw()
                win.flip()
            else:
                target.draw()
                pen.draw()
                train_aim_dot.draw()
                poke_dot.draw()
                win.flip() 
                

            # check for quit (the Esc key)
            if keyboard.state.has_key('escape'):
                win.close()
                core.quit()
            
        #in case there was no poke in the previous trial, add empty values to RT and poke position lists
        if poked: 
            poke_distance = dist(train_aim_dot.pos, poke_dot.pos)
            poke_pos      = poke_dot.pos
            poke_pos      = poke_pos.tolist()
        else:
            poke_distance    = float('nan')
            poke_time_train  = float('nan')
            poke_pos         = float('nan') 
      
        
        #increment after each trial
        if not released_too_soon:
            current_train_trial += 1
            #add training data to lists
            pokeDistancesTrain.append(poke_distance)
            train_aim_dot_position = train_aim_dot.pos.tolist()
            aimDotPositions.append(train_aim_dot_position)
            pokePositionsTrain.append(poke_pos)
            trialNumsTrain.append(current_train_trial)
            pokeTimesTrain.append(poke_time_train)
            releaseTimesTrain.append(release_time_train)

            #create dictionary from training data lists and dataframe from this dictionary
            train1_data_dict  = {'aimDotPos': aimDotPositions, 'pokePos': pokePositionsTrain, 'trial': trialNumsTrain, 'releaseTime': releaseTimesTrain, 'pokeTime': pokeTimesTrain, 'pokeDist': pokeDistancesTrain}
            train1_data_frame = pd.DataFrame(dict([(k,Series(v)) for k,v in train1_data_dict.iteritems()]))
            train1_data_frame.to_csv(_thisDir + '/data/%s_temp_reach_train1_output.csv' % (expInfo['participant']), mode = 'w')


    current_train_block += 1



#start out with zero total points
total_points = 0

current_train2_block = 0

#train2 loop
while current_train2_block < TOTAL_TRAIN2_BLOCK:
    # start block with zero block points
    block_points = 0


    train2_instr_timer.reset()
    train2_instr_timer.add(TRAIN2_INSTR_DUR)
    if current_train2_block == 0:
        train2_instr_text = 'PLEASE COME NOTIFY THE RA THAT YOU ARE FINISHED WITH PART 1.' 
    elif current_train2_block <= TOTAL_TRAIN2_BLOCK - 3: 
        train2_instr_text = 'You are about to begin block ' + str(current_train2_block + 1) + ' of ' + str(TOTAL_TRAIN2_BLOCK) + ' of part 2 of the training. Please press Z to continue.'
    else: #i.e., if they're in the last two blocks
        train2_instr_text = 'You are about to begin block ' + str(current_train2_block + 1) + ' of ' + str(TOTAL_TRAIN2_BLOCK) + ' of part 2 of the training. During the final two blocks, you will occasionally play for higher stakes. During these trials, the window will have a gold border, and the penalty values and reward values will be doubled. Please press Z to continue' 
    train2_instructions.setText(train2_instr_text)

    while not keyboard.state.has_key('z'):
        train2_instructions.draw()
        win.flip()
 
        # check for quit (the Esc key)
        if keyboard.state.has_key('escape'):
            win.close()
            core.quit()

    #TRIAL LOOP
    #duration of train2
    current_train2_trial = int(TRAIN2_TRIALS_PER_BLOCK * current_train2_block)
    total_train2_trial   = current_train2_trial + TRAIN2_TRIALS_PER_BLOCK
    
    while current_train2_trial < total_train2_trial:
        #------Prepare to start Routine "train2"-------
        poke_time_train2 = 0
        t = 0
        train2Clock.reset()  # clock 
        
        #Get and set positioning 
        sep_dis    = train2Info.iloc[current_train2_trial]['sepDistance']
        targ_pos_x = train2Info.iloc[current_train2_trial]['x_position']
        targ_pos_y = train2Info.iloc[current_train2_trial]['y_position']
        target.pos  = [targ_pos_x, targ_pos_y]
        penalty.pos = [targ_pos_x + sep_dis, targ_pos_y] # no magic numbers


        #Get and set reward/penalty values
        rew_value = train2Info.iloc[current_train2_trial]['reward']
        pen_value = train2Info.iloc[current_train2_trial]['penalty']

        if pen_value == 0:
            penalty.lineColor       = YELLOW
            penalty.fillColor       = GREY
            penalty_color.lineColor = YELLOW
            penalty_color.fillColor = YELLOW
        elif pen_value == -1:
            penalty.lineColor       = LBLUE
            penalty.fillColor       = LBLUE
            penalty_color.lineColor = LBLUE
            penalty_color.fillColor = LBLUE
        elif pen_value == -5: 
            penalty.lineColor       = RED
            penalty.fillColor       = RED 
            penalty_color.lineColor = RED
            penalty_color.fillColor = RED
        elif pen_value == -3:
            penalty.lineColor       = PURPLE
            penalty.fillColor       = PURPLE
            penalty_color.lineColor = PURPLE
            penalty_color.fillColor = PURPLE
            high_stakes_border.lineColor = PURPLE
        elif pen_value == -15:
            penalty.lineColor       = ORANGE
            penalty.fillColor       = ORANGE
            penalty_color.lineColor = ORANGE
            penalty_color.fillColor = ORANGE
            high_stakes_border.lineColor = ORANGE
        
        poked    = False
        poke_now = False
        released_too_soon = False
        space_pressed     = False

        #change fixation color to WHITE
        cross_vertical.fillColor   = WHITE
        cross_horizontal.fillColor = WHITE
        cross_vertical.lineColor   = WHITE
        cross_horizontal.lineColor = WHITE

        #reward/penalty display
        value_display_timer.reset()
        value_display_timer.add(VALUE_DISPLAY_DUR)
        penalty_value.setText('$ '+ str(pen_value))
        reward_value.setText('$' + str(rew_value))
        
        while value_display_timer.getTime() > 0:

            # check for quit (the Esc key)
            if keyboard.state.has_key('escape'):
                win.close()
                core.quit()

            if pen_value == -15 or pen_value == -3:
                high_stakes_border.draw()
                high_stakes_note_R.draw()
                high_stakes_note_L.draw()


            penalty_value.draw()
            reward_value.draw()
            penalty_color.draw()
            reward_color.draw()
           
            win.flip()
        
        space_bar_timer.reset()
        space_bar_timer.add(6)
        #FIXATION CROSS
        #Ends when spacebar is pressed
        while not keyboard.state.has_key(' '):
            # train2_text.setText(str(keyboard.state))
            # train2_text.draw()

              # check for quit (the Esc key)
            if keyboard.state.has_key('escape'):
                win.close()
                core.quit()
         
            if space_bar_timer.getTime > 0:
                press_text.draw()
            else:
                please_press.draw()

            cross_vertical.draw()
            cross_horizontal.draw()
            press_text.draw()
            

            win.flip()

            
        press_time = train2Clock.getTime()

        blue_fix_timer.reset()
        blue_fix_timer.add(random.uniform(.4, .6))
        #make fixation cross blue
        cross_vertical.fillColor   = BLUE
        cross_horizontal.fillColor = BLUE
        cross_vertical.lineColor   = BLUE
        cross_horizontal.lineColor = BLUE
        #BLUE FIXATION CROSS
        #Ends after time period between 400 and 600 ms
        while blue_fix_timer.getTime() > 0:

            
            cross_vertical.draw()
            cross_horizontal.draw()
            bound_box.draw()
            
            win.flip()
            

            if not keyboard.state.has_key(' '):
                released_too_soon = True

            # check for quit (the Esc key)
            if keyboard.state.has_key('escape'):
                win.close()
                core.quit()

        
        routineTimer.reset()
        routineTimer.add(TRIAL_TIME_LIMIT)
        continueRoutine = True
        train2Clock.reset()

        #initialize mouse position to bottom center of screen
        mouse.setPos(BOT_CENTER)


        target.opacity  = 1
        penalty.opacity = 1

        release_marked = False
        too_slow = False

        release_time_train2 = 0

        while not poked and not released_too_soon:
            
            # get current time
            t = TRIAL_TIME_LIMIT - routineTimer.getTime()
            frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
           
            if routineTimer.getTime() <= 0:
                too_slow = True
            # check for quit (the Esc key)
            if keyboard.state.has_key('escape'):
                win.close()
                core.quit()
        
            # #true iff poke occured at this frame within boundaries of 'background', i.e., anywhere in window AND the spacebar was released
            # poke_now =

            


            released = not keyboard.state.has_key(' ')

            if released and not release_marked:
                release_time_train2   = t
                release_marked        = True
            
            #too soon condition, 100ms
            if release_time_train2 < .1 and release_time_train2 != 0:
                released_too_soon = True
                break

            #if poke occurred at this frame, then prepare to draw poke_dot at that the poke location next frame and store time elapsed since trial start
            if mouse.isPressedIn(background, buttons=[0]) and  not keyboard.state.has_key(' '): 
                poked        = True
                poke_pos     = mouse.getPos()
                poke_dot.pos = poke_pos
                poke_time_train2  = t
                break

            background.draw()
            penalty.draw()
            target.draw()
            win.flip()

            

    

        #in case there was no poke in the previous trial, add empty values to RT, poke position, and distance lists
        if poked:
            penalty_distance = dist(penalty.pos, poke_dot.pos)
            target_distance  = dist(target.pos, poke_dot.pos)
            poke_pos         = poke_pos.tolist()
        else: 
            poke_time_train2 = float('nan')
            poke_pos         = float('nan')
            penalty_distance = float('nan')
            target_distance  = float('nan')

        poke_type = getPokeType(penalty_distance, target_distance, released_too_soon)
        points    = getTrialScore(poke_type, rew_value, pen_value)
        accuracy  = getTrialAcc(poke_type)

        if too_slow:
            points = TOO_SLOW_PEN
        
        point_text = '$ ' + str(points) 

        if points == 3:
            trial_point_feed.setColor(GOLD, colorSpace = 'rgb')
        elif points > 0:
            trial_point_feed.setColor(GREEN, colorSpace = 'rgb')
        elif points < 0:          
            trial_point_feed.setColor(RED, colorSpace = 'rgb')
        else:
            trial_point_feed.setColor(YELLOW, colorSpace = 'rgb')



        continueRoutine = True
        routineTimer.reset()  # clock 
        routineTimer.add(TRAIN2_FEED_DUR)

        frameCount = 0
        

        trial_point_feed.setText(point_text)

        #FEEDBACK LOOP
        while routineTimer.getTime() > 0:

            # get current time
            t =  TRAIN2_FEED_DUR - routineTimer.getTime()
            
            frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
            
            # change opacity of the circles every 10 frames, to give the appearance of flashing, based on where poke occurred
            frameCount = int(frameN/10)
            
        
            if released_too_soon:
                too_soon_text.draw()
                win.flip()
            elif too_slow:
                no_poke_text.draw()
                trial_point_feed.draw()
                # refresh the screen
                win.flip()
            else:
                opacities =  getOpacities(poke_type, frameCount)
                target.opacity  = opacities[0]
                penalty.opacity = opacities[1]
                penalty.draw() 
                target.draw() 
                poke_dot.draw()
                trial_point_feed.draw()
                # refresh the screen
                win.flip()
        

            # check for quit (the Esc key)
            if keyboard.state.has_key('escape'):
                win.close()
                core.quit()

        if not released_too_soon:
            #convert np.ndarray positions to lists
            target_position  = target.pos.tolist()
            penalty_position = penalty.pos.tolist()
            #add data to lists
            penaltyDistances.append(penalty_distance)
            targetDistances.append(target_distance)
            scores.append(points)
            pokePositionstrain2.append(poke_pos)
            penaltyPositions.append(penalty_position)
            targetPositions.append(target_position)
            pokeTimestrain2.append(poke_time_train2) #time distance from beginning of trial to poke
            #totalScores.append(total_points)
            trial_in_block = total_train2_trial - current_train2_trial - 1 #trial number within current block (0-17)
            trialNumber.append(trial_in_block)
            targetSizes.append(target.size[0]) #diameter
            blockNumber.append(current_train2_block)
            rewardValue.append(rew_value)
            penaltyValue.append(pen_value)
            sepDistances.append(sep_dis)
            timeLimits.append(TRIAL_TIME_LIMIT)
            accuracies.append(accuracy)
            pokeTypes.append(poke_type)
            releaseTimesTrain2.append(release_time_train2)

                        #create dictionary from train2 data lists and dataframe from this dictionary
            train2_data_dict  = {'score': scores, 'pokePos': pokePositionstrain2, 'penaltyPos': penaltyPositions, 'targDist': targetDistances,
                               'penalDist': penaltyDistances, 'targetPos': targetPositions, 'pokeTime': pokeTimestrain2, 'releaseTime': releaseTimesTrain2, 'trial': trialNumber, 'block': blockNumber, 'targSize': targetSizes, 
                               'sepDist' : sepDistances, 'rewardVal': rewardValue, 'penaltyVal': penaltyValue, 'accuracy': accuracies, 'pokeType': pokeTypes}
            train2_data_frame = pd.DataFrame(dict([(k,Series(v)) for k,v in train2_data_dict.iteritems()]))
            train2_data_frame.to_csv(_thisDir + '/data/%s_temp_reach_train2_output.csv' % (expInfo['participant']), mode = 'w')


            #increment block points after each successful trial
            block_points += points

            #increment to next trial
            current_train2_trial += 1

        else:#if they did release to soon, place this trial at the end of the block to be tried again
            if current_train2_trial == 0:
                train2Info = pd.concat([train2Info.ix[1:total_train2_trial - 1], train2Info.ix[[0]], train2Info.ix[total_train2_trial:]])
                train2Info = train2Info.reset_index(drop=True)
            elif total_train2_trial < len(train2Info.index) and not (current_train2_trial + 1) % TRAIN2_TRIALS_PER_BLOCK == 0:
                train2Info = pd.concat([train2Info.ix[0:current_train2_trial - 1], 
                            train2Info.ix[current_train2_trial + 1:total_train2_trial - 1], train2Info.ix[[current_train2_trial]], train2Info.ix[total_train2_trial:]])
                train2Info = train2Info.reset_index(drop=True)
            elif not (current_train2_trial + 1) % TRAIN2_TRIALS_PER_BLOCK == 0:
                train2Info = pd.concat([train2Info.ix[0:current_train2_trial - 1], train2Info.ix[current_train2_trial + 1:], train2Info.ix[[current_train2_trial]]])
                train2Info = train2Info.reset_index(drop=True)
    
    total_points += block_points

    current_train2_block += 1

    #timed feedback slide after each block
    #point_feed_timer.reset()
    #point_feed_timer.add(POINT_FEED_DUR)
    #block_point_feed.setText('Points in previous block: ' + str(block_points))
    #total_point_feed.setText('Total points: ' + str(total_points))

    #while point_feed_timer.getTime() > 0:

    #    encourage_feed.draw()
    #    block_point_feed.draw()
    #    total_point_feed.draw()
    #    win.flip()
        
train2_data_frame.to_csv(_thisDir + '/data/%s_reach_train2_output.csv' % (expInfo['participant']), mode = 'w')
train1_data_frame.to_csv(_thisDir + '/data/%s_reach_train1_output.csv' % (expInfo['participant']), mode = 'w')

   
win.close()
core.quit()
