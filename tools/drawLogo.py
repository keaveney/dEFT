# standard imports
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# displays logos inline within the notebook;
# remove if using a python interpreter instead
#%matplotlib inline

# logomaker import
import logomaker


# make Figure and Axes objects
fig, ax = plt.subplots(1,1,figsize=[4,2])

# load logo matrix
logo_df = logomaker.get_example_matrix('logomaker_logo_matrix',
                                       print_description=True)


print(logo_df)

#data = [['Alex',10],['Bob',12],['Clarke',13]]

data = [[0.8,  0.0, 0.0, -0.5,  0.0,  0.0,  0.0  0.0],[ 0.0  0.6  0.0  0.0 -0.5  0.0  0.0  0.0],[ 0.0  0.0  0.6  0.0  0.0 -0.5  0.0  0.0], [0.0  0.6  0.0  0.0  0.0  0.0 -0.5  0.0],[0.0  0.0  0.0  0.0  0.0  0.0  0.0 -0.5]]
df = pd.DataFrame(data,columns=['d','E','F','T'])
print(df)

#logo_df =
#       L    O    G    m    a    k    e    r
#pos
#0    0.8  0.0  0.0 -0.5  0.0  0.0  0.0  0.0
#1    0.0  0.6  0.0  0.0 -0.5  0.0  0.0  0.0
#2    0.0  0.0  0.6  0.0  0.0 -0.5  0.0  0.0
#3    0.0  0.6  0.0  0.0  0.0  0.0 -0.5  0.0
#4    0.0  0.0  0.0  0.0  0.0  0.0  0.0 -0.5

# create color scheme
color_scheme = {
    'd' : [0, .5, 0],
    'E' : [1, 0, 0],
    'F' : [1, .65, 0],
    'T' : [1, .65, 0],
    'maker': 'gray'
}

# create Logo object
logo_logo = logomaker.Logo(logo_df,
                           ax=ax,
                           color_scheme=color_scheme,
                           baseline_width=0,
                           font_name='Arial',
                           show_spines=False,
                           vsep=.005,
                           width=.95)

# color the 'O' at the end of the logo a different color
logo_logo.style_single_glyph(c='O', p=3, color=[0, 0, 1])

# change the font of 'maker' and flip characters upright.
logo_logo.style_glyphs_below(font_name='OCR A Std', flip=False, width=1.0)

# remove tick marks
ax.set_xticks([])
ax.set_yticks([])

# tighten layout
logo_logo.fig.tight_layout()

fig.savefig("logo.png")
