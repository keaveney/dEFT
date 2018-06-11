from asciimatics.effects import Cycle, Stars
from asciimatics.renderers import FigletText
from asciimatics.scene import Scene
from asciimatics.screen import Screen

def deft_splash(screen):
    effects = [
        Cycle(
            screen,
            FigletText("dEFT", font='big'),
            int(screen.height / 2 - 8)),
        Cycle(
            screen,
            FigletText("Differential Effective Field Theory Fitter", font='small'),
            int(screen.height / 2 + 3)),
        Stars(screen, 900)
    ]
    screen.play([Scene(effects, 900)])

Screen.wrapper(deft_splash)
