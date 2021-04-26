import random
import pandas as pd

NUM_OF_BLOCKS = 8
COPIES_PER_BLOCK = 3
TIME_LIMIT = 1 #seconds
SMALL_REWARD = 1
HUGE_REWARD = 3 * SMALL_REWARD
SMALL_PENALTY = 0
MED_PENALTY = -1
LARGE_PENALTY = -5
HUGE_PENALTY = 3 * MED_PENALTY
HUGE_PENALTY2 = 3 * LARGE_PENALTY
RADIUS = 32  #pixels
SMALL_SEP_LEFT = -RADIUS
MED_SEP_LEFT = -1.375 * RADIUS
#LARGE_SEP_LEFT = -1.5 * RADIUS
SMALL_SEP_RIGHT = RADIUS
MED_SEP_RIGHT = 1.375 * RADIUS
#LARGE_SEP_RIGHT = 1.5 * RADIUS
Loc1 = [-100, -100, False, False, False, False, False] #[xPos, yPos, loc1 used recently for pen = 0, used recently for pen = -100, used recently for pen = -500...]
Loc2 = [-100, 0, False, False, False, False, False]
Loc3 = [-100, 100, False, False, False, False, False]
Loc4 = [0, -100, False, False, False, False, False]
Loc5 = [0, 0, False, False, False, False, False]
Loc6 = [0, 100, False, False, False, False, False]
Loc7 = [100, -100, False, False, False, False, False]
Loc8 = [100, 0, False, False, False, False, False]
Loc9 = [100, 100, False, False, False, False, False]
CONFIG = {"blockNum":[], "reward":[], "penalty":[], "sepDistance":[], "timeLimit":[], "x_position":[], "y_position":[]} #dictionary holding all possible configurations

#This function searches through list of 9 possible locations and assigns a location randomly to configuration in the dataframe. 
#Only locations which have not been used recently for the specific configuration will be assigned
def insertLocation(daData, locList, pen, row):
    found = False
    if pen == SMALL_PENALTY:
        n = 0
        while not found:
            locNum = random.randint(0, 8)
            if not locList[locNum][2]:
                found = True
                daData.loc[(row, "x_position")] = locList[locNum][0]
                daData.loc[(row, "y_position")] = locList[locNum][1]
                locList[locNum][2] = True
            n += 1
            if n == 8:
                resetLocValues(locList, 2)
    elif pen == MED_PENALTY:
        n = 0
        while not found:
            locNum = random.randint(0, 8)
            if not locList[locNum][3]:
                found = True
                daData.loc[(row, "x_position")] = locList[locNum][0]
                daData.loc[(row, "y_position")] = locList[locNum][1]
                locList[locNum][3] = True
            n += 1
            if n == 8:
                resetLocValues(locList, 3)
    elif pen == LARGE_PENALTY:
        n = 0
        while not found:
            locNum = random.randint(0, 8)
            if not locList[locNum][4]:
                found = True
                daData.loc[(row, "x_position")] = locList[locNum][0]
                daData.loc[(row, "y_position")] = locList[locNum][1]
                locList[locNum][4] = True
            n += 1
            if n == 8:
                resetLocValues(locList, 4)
    elif pen == HUGE_PENALTY:
        n = 0
        while not found:
            locNum = random.randint(0, 8)
            if not locList[locNum][5]:
                found = True
                daData.loc[(row, "x_position")] = locList[locNum][0]
                daData.loc[(row, "y_position")] = locList[locNum][1]
                locList[locNum][5] = True
            n += 1
            if n == 8:
                resetLocValues(locList, 5)
    elif pen == HUGE_PENALTY2:
        n = 0
        while not found:
            locNum = random.randint(0, 8)
            if not locList[locNum][6]:
                found = True
                daData.loc[(row, "x_position")] = locList[locNum][0]
                daData.loc[(row, "y_position")] = locList[locNum][1]
                locList[locNum][6] = True
            n += 1
            if n == 8:
                resetLocValues(locList, 6)

#this reset all the "has been used recently" values in the Locations list to false 
def resetLocValues(List, pos):
    for x in range(0, len(List)):
        List[x][pos] = False

#This function creates all possible configurations of penalty and sep distances and stores the in config list
def createConfigurations(configList):
    for x in range (0, 3): #possible penalty values
        for y in range (0, 4): #4 possible seperation distances
            
            #append blockNum
            configList["blockNum"].append(0)

            #append reward
            configList["reward"].append(SMALL_REWARD)

            #append penalty
            if x == 0:
                configList["penalty"].append(SMALL_PENALTY)
            elif x == 1:
                configList["penalty"].append(MED_PENALTY)
            elif x == 2:
                configList["penalty"].append(LARGE_PENALTY)


            #append separation distance
            if y == 0:
                configList["sepDistance"].append(SMALL_SEP_LEFT)
            elif y == 1:
                configList["sepDistance"].append(SMALL_SEP_RIGHT)
            elif y == 2:
                configList["sepDistance"].append(MED_SEP_LEFT)
            else:
                configList["sepDistance"].append(MED_SEP_RIGHT)
            

            #append timeLimit
            configList["timeLimit"].append(TIME_LIMIT)

            #append filler positions, will be updated later 
            configList["x_position"].append(0)
            configList["y_position"].append(0)

def createLocList(locList):
    locList.append(Loc1)
    locList.append(Loc2)
    locList.append(Loc3)
    locList.append(Loc4)
    locList.append(Loc5)
    locList.append(Loc6)
    locList.append(Loc7)
    locList.append(Loc8)
    locList.append(Loc9)

def addSpecialRewards(df, numRewards):
    for x in range(0, 2):
        for y in range(0, numRewards):
            row = len(df.index) - random.randint(1, COPIES_PER_BLOCK * 12)
            if x == 0 and y % 2:
                df.loc[(row, "reward")] = HUGE_REWARD
                df.loc[(row, "penalty")] = HUGE_PENALTY
            elif x == 0:
                df.loc[(row, "reward")] = HUGE_REWARD
                df.loc[(row, "penalty")] = HUGE_PENALTY2
            elif x == 1 and y % 2:
                row -= COPIES_PER_BLOCK * 12
                df.loc[(row, "reward")] = HUGE_REWARD
                df.loc[(row, "penalty")] = HUGE_PENALTY
            else:
                row -= COPIES_PER_BLOCK * 12
                df.loc[(row, "reward")] = HUGE_REWARD
                df.loc[(row, "penalty")] = HUGE_PENALTY2


def main():
    #create all possible configurations
    createConfigurations(CONFIG)

    #load all dataFrames to be concatenated into list. Shuffle each dataframe prior to loading
    toConcat = []
    for y in range(0, NUM_OF_BLOCKS):
        for z in range(0, len(CONFIG["blockNum"])):
            CONFIG["blockNum"][z] = y 
        data = pd.DataFrame(CONFIG)
        for x in range(0, COPIES_PER_BLOCK):
            data = data.sample(frac = 1).reset_index(drop=True)
            toConcat.append(data)

    #concatenate list of dataframes into final dataframe and reset indexing
    finalData = pd.concat(toConcat)
    finalData = finalData.reset_index(drop=True)

    #create possible positions from which to choose locations randomly
    Locations = []
    createLocList(Locations)
    
    numSpecialRewards = 12
    addSpecialRewards(finalData, numSpecialRewards)

    #insert locations into every row
    for x in range(0, len(finalData.index)):
        if finalData.iloc[x]["penalty"] == LARGE_PENALTY:
            insertLocation(finalData, Locations, LARGE_PENALTY, x)
        elif finalData.iloc[x]["penalty"] == MED_PENALTY:
            insertLocation(finalData, Locations, MED_PENALTY, x)
        elif finalData.iloc[x]["penalty"] == SMALL_PENALTY:
            insertLocation(finalData, Locations, SMALL_PENALTY, x)
        elif finalData.iloc[x]["penalty"] == HUGE_PENALTY:
            insertLocation(finalData, Locations, HUGE_PENALTY, x)
        elif finalData.iloc[x]["penalty"] == HUGE_PENALTY2:
            insertLocation(finalData, Locations, HUGE_PENALTY2, x)

    

    #print(finalData)

    finalData.to_csv('trainTrialsPt2.csv')

    #print(pd.read_csv('test.csv'))

main()