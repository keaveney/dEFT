from asciimatics.effects import Cycle, Stars
from asciimatics.renderers import FigletText
from asciimatics.scene import Scene
from time import sleep

def deft_splash(screen):
    effects = [
               Cycle(
                     screen,
                     FigletText("dEFT", font='big'),
                     int(screen.height / 2 - 8)),
               Cycle(
                     screen,
                     FigletText("A differential Effective Field Theory tool", font='small'),
                     int(screen.height / 2 + 3)),
               Cycle(
                     screen,
                     FigletText("press q to run job!", font='small'),
                     int(screen.height / 1.5)),
               Stars(screen, 900)
               ]
    screen.play([Scene(effects, 900)])
    ev = screen.get_key()
    if ev in (ord('Q'), ord('q')):
            return
