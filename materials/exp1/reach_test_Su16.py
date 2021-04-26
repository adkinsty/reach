
from reach_constants_Su16 import *
from reach_functions_Su16 import *
from reach_components_Su16 import *

# Ensure that relative paths start from the same directory as this script
_thisDir = os.path.dirname(os.path.abspath(__file__)).decode(sys.getfilesystemencoding())
os.chdir(_thisDir)

#start out with zero total points
total_points = 0

current_test_block = 0

#mouse.setVisible(0)
#Test loop
while current_test_block < TOTAL_TEST_BLOCK:
    # start block with zero block points
    block_points = 0


    test_instr_timer.reset()
    test_instr_timer.add(TEST_INSTR_DUR)
    test_instr_text = 'You are about to begin block ' + str(current_test_block + 1) + ' of ' + str(TOTAL_TEST_BLOCK) + '.\nTry to earn as many points as you can by poking within the GREEN target circle with your finger. Please press Z to continue.'
    test_instructions.setText(test_instr_text)

    while not keyboard.state.has_key('z'):
        test_instructions.draw()
        win.flip() 
 
        # check for quit (the Esc key)
        if keyboard.state.has_key('escape'):
            win.close()
            core.quit()

    #TRIAL LOOP
    #duration of test
    current_test_trial = int(TEST_TRIALS_PER_BLOCK * current_test_block)
    total_test_trial = current_test_trial + TEST_TRIALS_PER_BLOCK
    
    while current_test_trial < total_test_trial:
        #------Prepare to start Routine "Test"-------
        poke_time_test = 0
        t = 0
        TestClock.reset()  # clock 
        frameN = -1
		
        #Get and set positioning 
        sep_dis    = testInfo.iloc[current_test_trial]['sepDistance']
        targ_pos_x = testInfo.iloc[current_test_trial]['x_position']
        targ_pos_y = testInfo.iloc[current_test_trial]['y_position']
        target.pos  = [targ_pos_x, targ_pos_y]
        penalty.pos = [targ_pos_x + sep_dis, targ_pos_y] # no magic numbers


        #Get and set reward/penalty values
        rew_value = testInfo.iloc[current_test_trial]['reward']
        pen_value = testInfo.iloc[current_test_trial]['penalty']

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
        
        #elif pen_value == -2:
        #elif pen_value == -10:
        
        poked    = False
        poke_now = False
        released_too_soon = False
       
        #-------Start Routine "Test"-------

        space_pressed = False

        #change fixation color to WHITE
        cross_vertical.fillColor   = WHITE
        cross_horizontal.fillColor = WHITE
        cross_vertical.lineColor   = WHITE
        cross_horizontal.lineColor = WHITE

        #reward/penalty display
        value_display_timer.reset()
        value_display_timer.add(VALUE_DISPLAY_DUR)
        penalty_value.setText('$ ' + str(pen_value))
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


            
        #FIXATION CROSS
        #Ends when spacebar is pressed
        space_bar_timer.reset()
        space_bar_timer.add(6)
        while not keyboard.state.has_key(' '):
            # test_text.setText(str(keyboard.state))
            # test_text.draw()

              # check for quit (the Esc key)
            if keyboard.state.has_key('escape'):
                win.close()
                core.quit()

            if space_bar_timer.getTime() > 0:
                press_text.draw()
            else:    
                please_press.draw()
            cross_vertical.draw()
            cross_horizontal.draw()
            
            win.flip()
            
        press_time = TestClock.getTime()

        blue_fix_timer.reset()
        blue_fix_timer.add(random.uniform(.4, .6))
        
        #BLUE FIXATION CROSS
        #Ends after time period between 400 and 600 ms
        while blue_fix_timer.getTime() > 0:

            #make fixation cross blue
            cross_vertical.fillColor   = BLUE
            cross_horizontal.fillColor = BLUE
            cross_vertical.lineColor   = BLUE
            cross_horizontal.lineColor = BLUE

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
        TestClock.reset()

        #initialize mouse position to bottom center of screen
        mouse.setPos(BOT_CENTER)
        
        target.opacity  = 1
        penalty.opacity = 1

        release_marked = False
        too_slow = False

        release_time_test = 0

        while not poked and not released_too_soon:
            # get current time
            t = TRIAL_TIME_LIMIT - routineTimer.getTime()

            if routineTimer.getTime() <= 0:
                too_slow = True


            frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
           
            
            
            
            # check for quit (the Esc key)
            if keyboard.state.has_key('escape'):
                win.close()
                core.quit()
        
            #true iff poke occured at this frame within boundaries of 'background', i.e., anywhere in window AND the spacebar was released
            # poke_now = mouse.isPressedIn(background, buttons=[0]) and  not keyboard.state.has_key(' ')



            released = not keyboard.state.has_key(' ')

            if released and not release_marked:
                release_time_test   = t
                release_marked       = True
            
            #too soon condition, 100ms
            if release_time_test < .1 and release_time_test != 0:
                released_too_soon = True
                break

            #if poke occurred at this frame, then prepare to draw poke_dot at that the poke location next frame and store time elapsed since trial start
            if mouse.isPressedIn(background, buttons=[0]) and  not keyboard.state.has_key(' '): 
                poked        = True
                poke_pos     = mouse.getPos()
                poke_dot.pos = poke_pos
                poke_time_test  = t
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
            poke_time_test   = float('nan')
            poke_pos         = float('nan')
            penalty_distance = float('nan')
            target_distance  = float('nan')

        poke_type = getPokeType(penalty_distance, target_distance, released_too_soon)
        points    = getTrialScore(poke_type, rew_value, pen_value)
        accuracy  = getTrialAcc(poke_type)

        if too_slow:
        	points = TOO_SLOW_PEN
        

        point_text = '$ ' + str(points) 

        if points == 1:
            trial_point_feed.setColor(GREEN, colorSpace = 'rgb')
        elif points < 0:
            trial_point_feed.setColor(RED, colorSpace = 'rgb')
        elif points == 3:
            trial_point_feed.setColor(GOLD, colorSpace = 'rgb')
        else:
            trial_point_feed.setColor(YELLOW, colorSpace = 'rgb')



        continueRoutine = True
        routineTimer.reset()  # clock 
        routineTimer.add(TEST_FEED_DUR)

        frameCount = 0
        

        trial_point_feed.setText(point_text)

        #FEEDBACK LOOP
        while routineTimer.getTime() > 0:

            # get current time
            t =  TEST_FEED_DUR - routineTimer.getTime()
            
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
            too_slow_bools.append(too_slow)
            penaltyDistances.append(penalty_distance)
            targetDistances.append(target_distance)
            scores.append(points)
            pokePositionsTest.append(poke_pos)
            penaltyPositions.append(penalty_position)
            targetPositions.append(target_position)
            pokeTimesTest.append(poke_time_test) #time distance from beginning of trial to poke
            releaseTimesTest.append(release_time_test)
            
            #totalScores.append(total_points)
            trial_in_block = total_test_trial - current_test_trial - 1 #trial number within current block (0-17)
            trialNumber.append(trial_in_block)
            targetSizes.append(target.size[0]) #diameter
            blockNumber.append(current_test_block)
            rewardValue.append(rew_value)
            penaltyValue.append(pen_value)
            sepDistances.append(sep_dis)
            timeLimits.append(TRIAL_TIME_LIMIT)
    
            accuracies.append(accuracy)
            pokeTypes.append(poke_type)

                        #create dictionary from test data lists and dataframe from this dictionary
            test_data_dict  = {'tooSlow':too_slow_bools, 'score': scores, 'pokePos': pokePositionsTest, 'penaltyPos': penaltyPositions, 'targDist': targetDistances,
                               'penalDist': penaltyDistances, 'targetPos': targetPositions, 'pokeTime': pokeTimesTest, 'releaseTime': releaseTimesTest, 'trial': trialNumber, 'block': blockNumber, 'targSize': targetSizes, 
                               'sepDist' : sepDistances, 'rewardVal': rewardValue, 'penaltyVal': penaltyValue, 'accuracy': accuracies, 'pokeType': pokeTypes}
            test_data_frame = pd.DataFrame(dict([(k,Series(v)) for k,v in test_data_dict.iteritems()]))
            test_data_frame.to_csv(_thisDir + '/data/%s_temp_reach_test_output.csv' % (expInfo['participant']), mode = 'w')


            #increment block points after each successful trial
            block_points += points

            #increment to next trial
            current_test_trial += 1

        else:#if they did release to soon, place this trial at the end of the block to be tried again
            if current_test_trial == 0:
                testInfo = pd.concat([testInfo.ix[1:total_test_trial - 1], testInfo.ix[[0]], testInfo.ix[total_test_trial:]])
                testInfo = testInfo.reset_index(drop=True)
            elif total_test_trial < len(testInfo.index) and not (current_test_trial + 1) % TEST_TRIALS_PER_BLOCK == 0:
                testInfo = pd.concat([testInfo.ix[0:current_test_trial - 1], 
                            testInfo.ix[current_test_trial + 1:total_test_trial - 1], testInfo.ix[[current_test_trial]], testInfo.ix[total_test_trial:]])
                testInfo = testInfo.reset_index(drop=True)
            elif not (current_test_trial + 1) % TEST_TRIALS_PER_BLOCK == 0:
                testInfo = pd.concat([testInfo.ix[0:current_test_trial - 1], testInfo.ix[current_test_trial + 1:], testInfo.ix[[current_test_trial]]])
                testInfo = testInfo.reset_index(drop=True)
    
    total_points += block_points

    current_test_block += 1

        
test_data_frame.to_csv(_thisDir + '/data/%s_reach_test_output.csv' % (expInfo['participant']), mode = 'w')


bonuses        = []
selectedTrials = []
for b in range(TOTAL_TEST_BLOCK):
    selected_trial = TEST_TRIALS_PER_BLOCK * b + randint(0,TEST_TRIALS_PER_BLOCK)
    bonus = test_data_frame.iloc[selected_trial]['score']
    bonuses.append(bonus)
    selectedTrials.append(selected_trial)

total_bonus   = sum(bonuses)    
payment_dict  = {'trials': selectedTrials, 'bonuses': bonuses}  
payment_frame = test_data_frame = pd.DataFrame(dict([(k,Series(v)) for k,v in payment_dict.iteritems()]))

print payment_frame
print 'Your cash bonus is: $' + str(total_bonus)


win.close()
core.quit()