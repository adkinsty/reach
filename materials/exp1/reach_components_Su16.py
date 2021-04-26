import sys
ipy_stdout = sys.stdout
from psychopy import locale_setup, visual, core, data, event, logging, sound, gui
from psychopy.constants import *  # things like STARTED, FINISHED
from reach_constants_Su16 import *
import random
from pyglet.window import key, Window
from psychopy.iohub import launchHubServer
sys.stdout = ipy_stdout
# Start the ioHub process. The return variable is what is used
# during the experiment to control the iohub process itself,
# as well as any running iohub devices.
io=launchHubServer()

# Setup the Window
win = visual.Window(size=[900, 900], fullscr=False, screen=0, allowGUI=True, allowStencil=False,
    monitor=u'testMonitor', color=[0,0,0], colorSpace='rgb',
    blendMode='add', useFBO=True,
    units='pix')

# # store frame rate of monitor if we can measure it successfully
expInfo['frameRate']=win.getActualFrameRate()
if expInfo['frameRate']!=None:
    frameDur = 1.0/round(expInfo['frameRate'])
else:
    frameDur = 1.0/60.0 # couldn't get a reliable measure so guess

# By default, ioHub will create Keyboard and Mouse devices and
# start monitoring for any events from these devices only.
keyboard = io.devices.keyboard
mouse    = event.Mouse(visible=True, win=win)


#initialize transparent background shape; passing this shape is necessary for executing the mouse.isPressedIn() function
background = visual.Polygon(win=win, name='background',
	edges = 4, size=[2000, 2000],
	ori=0, pos=[0, 0],
	lineWidth=1, lineColor=[255, 255, 0], lineColorSpace='rgb',
	fillColor=[255, 255, 0], fillColorSpace='rgb',
	opacity=0,depth=-1.0, 
	interpolate=True)
high_stakes_border = visual.Polygon(win=win, name='border',
	edges = 4, size=[900, 900],
	ori=45, pos=[0, 0],
	lineWidth=20, lineColor= GOLD, lineColorSpace='rgb',
	fillColor=[0,0,0], fillColorSpace='rgb',
	opacity=1,depth=-1.0, 
	interpolate=True)
train_aim_dot = visual.Polygon(win=win, name='train_aim_dot',
	edges = 180, size=[10, 10],
	ori=0, pos=[0, 0],
	lineWidth=1, lineColor=[255, 255, 0], lineColorSpace='rgb',
	fillColor=[255, 255, 0], fillColorSpace='rgb',
	opacity=1,depth=-1.0, 
	interpolate=True)
poke_dot = visual.Polygon(win=win, name='poke_dot',
	edges = 180, size=[5, 5],
	ori=0, pos=[0, 0],
	lineWidth=1, lineColor=[-1, -1, 1], lineColorSpace='rgb',
	fillColor=[-1, -1, 1], fillColorSpace='rgb',
	opacity=1,depth=-1.0, 
	interpolate=True)
target = visual.Polygon(win=win, name='target',
	edges = 180,  size=[64, 64],
	ori=0, pos=[0, 0],
	lineWidth=2, lineColor=[-1, 1, -1], lineColorSpace='rgb',
	fillColor=None, fillColorSpace='rgb',
	opacity=1,depth=0.0, 
	interpolate=True, units = 'pix')
pen = visual.Polygon(win=win, name='pen',
	edges = 180,  size=[64, 64],
	ori=0, pos=[0, 0],
	lineWidth=2, lineColor=YELLOW, lineColorSpace='rgb',
	fillColor=None, fillColorSpace='rgb',
	opacity=1,depth=0.0, 
	interpolate=True, units = 'pix')
penalty_color = visual.Polygon(win=win, name='penalty_color', edges = 180, 
	size = [64, 64], ori = 0, pos = [0, -64], lineWidth = 5, 
	lineColor = [1, -1, -1], lineColorSpace = 'rgb', 
	fillColor = [1, -1, -1], fillColorSpace = 'rgb',
	opacity = 1, depth = -1, interpolate = True)
reward_color = visual.Polygon(win=win, name='reward_color', edges = 180, 
	size = [64, 64], ori = 180, pos = [0, 64], lineWidth = 5, 
	lineColor = [-1, 1, -1], lineColorSpace = 'rgb', 
	fillColor = [-1, 1, -1], fillColorSpace = 'rgb',
	opacity = 1, depth = -1, interpolate = True)
penalty = visual.Polygon(win=win, name='penalty',
	edges = 180,   size=[64, 64],
	ori=0, pos=[-32, 0],
	lineWidth=1, lineColor=[1, -1, -1], lineColorSpace='rgb',
	fillColor=None, fillColorSpace='rgb',
	opacity=1,depth=-1.0, 
	interpolate=True, units = 'pix')
cross_horizontal = visual.Polygon(win=win, name='horiz',
	edges = 2, size=[60, 60],
	ori=0, pos=[0, 0], 
	lineWidth=5, lineColor=[1,1,1], lineColorSpace='rgb',
	fillColor=[1, 1, 1], fillColorSpace='rgb',
	opacity=1,depth=-1.0, 
	interpolate=True)
cross_vertical = visual.Polygon(win=win, name='vertic',
	edges = 2, size=[60, 60],
	ori=90, pos=[0, 0],
	lineWidth=5, lineColor=[1,1,1], lineColorSpace='rgb',
	fillColor=[1, 1, 1], fillColorSpace='rgb',
	opacity=1,depth=-1.0, 
	interpolate=True)
bound_box = visual.Polygon(win=win, name='bound_box',
	edges = 4,  size=[500, 500],
	ori=45, pos=[0, 0],
	lineWidth=1, lineColor=[-1, -1, 1], lineColorSpace='rgb',
	fillColor=None, fillColorSpace='rgb',
	opacity=1,depth=-1.0, 
	interpolate=True, units = 'pix')


# bound_box.setSize([160, 160], units = 'pix')
# target.setSize([64, 64], units = 'pix')
# penalty.setSize([64, 64], units = 'pix')


#text objects
no_poke_text       = visual.TextStim(win, text = u"Too Slow!")
too_soon_text      = visual.TextStim(win, text = u"You Reached Too Soon \n\nWait For Target")
test_text    	   = visual.TextStim(win, text = u"")
train2_text        = visual.TextStim(win, text = u"")
press_text  	   = visual.TextStim(win, pos =(0,50), text='Hold down the SPACEBAR...')
penalty_value	   = visual.TextStim(win, pos =(0, -120), text=' ')
reward_value  	   = visual.TextStim(win, pos =(0, 120), text=' ')
block_point_feed   = visual.TextStim(win, pos = (0, -100), text = '')
total_point_feed   = visual.TextStim(win, pos = (0, -200), text = '')
encourage_feed     = visual.TextStim(win, pos = (0, 200), text =  '')
train_instructions = visual.TextStim(win, pos = (0, 0), text = '')
train_instructions.setSize(200, units='pixels')
test_instructions  = visual.TextStim(win, pos = (0, 50), text = 'Nice job! Check out your scores:')
train2_instructions= visual.TextStim(win, pos = (0, 50), text = 'Nice job! Check out your scores:')
trial_point_feed   = visual.TextStim(win, pos = (0, 200), text = '')
payment_feed       = visual.TextStim(win, pos = (0, 0), text = '')
high_stakes_note_R = visual.TextStim(win, pos = (500, 0), ori = 90, text = 'HIGH STAKES ROUND!', colorSpace = 'rgb', color = GOLD)
high_stakes_note_L = visual.TextStim(win, pos = (-500, 0), ori = -90, text = 'HIGH STAKES ROUND!', colorSpace = 'rgb', color = GOLD)  
please_press       = visual.TextStim(win, pos =(0,50), text='PLEASE HOLD DOWN THE SPACEBAR NOW; Remember, you can take breaks between each block.')


# Create some handy timers
globalClock          = core.Clock()  # to track the time since experiment started
routineTimer         = core.CountdownTimer()  # to track time remaining of each (non-slip) routine 
TestClock            = core.Clock()
train2Clock          = core.Clock()  
trainingClock        = core.Clock()
blue_fix_timer       = core.CountdownTimer()
aim_dot_timer        = core.CountdownTimer()
feedback_timer       = core.CountdownTimer()
value_display_timer  = core.CountdownTimer()
train_instr_timer    = core.CountdownTimer()
test_instr_timer     = core.CountdownTimer()
train2_instr_timer   = core.CountdownTimer()
feedback_timer       = core.CountdownTimer()
point_feed_timer     = core.CountdownTimer()
payment_feed_timer   = core.CountdownTimer()
space_bar_timer      = core.CountdownTimer()
