"""
Copyright (C) 2022, Manel Vila-Vidal
Contact details: m@vila-vidal.com / manel.vila-vidal@upf.edu
Date: 10 Mar 2022

-----------------------------------------------------------------------

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.

-----------------------------------------------------------------------

If you use the source code, please make sure to reference both the
package and the paper:

> Vila-Vidal, M. (2022). VisualContour v1.0,
https://github.com/compaes/Visual-Contour. Zenodo. (DOI)

> REFERENCE PAPER
"""


import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style='white',color_codes=True,font_scale=2)
from ipywidgets import interact, interactive, interactive_output, fixed, interact_manual
import ipywidgets as widgets

def cardinalCubicHermiteSpline(s,x,sinterp,t=0,ends=False):
    # t is the tension parameter
    [N,dim]=x.shape
    m=np.zeros(x.shape)
    k=0
    m[k]=(1-t)*(x[k+1]-x[k])/(s[k+1]-s[k])
    for k in range(1,N-1):
        m[k]=(1-t)*(x[k+1]-x[k-1])/(s[k+1]-s[k-1])
    k+=1
    m[k]=(1-t)*(x[k]-x[k-1])/(s[k]-s[k-1])
    
    Ninterp=sinterp.size
    xinterp=[]
    for i in range(Ninterp):
        ss=sinterp[i]
        k=np.where(s<=ss)[0][-1]
        if ((k==0)+(k==N-1)+(k==N-2))*(ends==False): continue
        if (k==N-1)*(ends==True):
            xinterp+=[x[-1]]
            break
        sk=(ss-s[k])/(s[k+1]-s[k])
        xinterp+=[(2*sk**3-3*sk**2+1)*x[k]  +  (sk**3-2*sk**2+sk)*m[k]  +  (-2*sk**3+3*sk**2)*x[k+1]  +  (sk**3-sk**2)*m[k+1]]
    
    return np.array(xinterp)

def handler(change):
    Randomize.value=False

def contourgeneratorinteractive(R,v,d,t):

    N=int(v/2)
    D=d*1.3
    
    global RAN
    if R: RAN=np.random.randn(3*Nmax)
    Randomize.observe(handler,'value')
    
    # point radius tolerance, points are placed at the ring ± a radius tolerance following a normal with std=er:
    er=0.2*D 
    # point angle tolerance, points are uniformly distributed along the circle ± an angle tolerance following a normal with std=ea:
    ea=0.2*2*np.pi/float(N) 

    # outer ring
    ro=(1.5+D/2.)+er*RAN[:N] # outer radiuses
    #ao=2*np.pi*np.random.rand()+np.linspace(0,2*np.pi*(N-1)/float(N),N)+ea*np.random.randn(N)
    ao=np.linspace(0,2*np.pi*(N-1)/float(N),N)+ea*RAN[N:2*N] # outer angles
    ao.sort()
    xo=ro*np.cos(ao)
    yo=ro*np.sin(ao)

    # inner ring
    ri=(1.5-d/2)+er*RAN[2*N:3*N] # inner radiuses
    aux=np.array(list(ao[1:])+[ao[0]+2*np.pi])
    ai=(aux+ao)/2 # inner angles = each point is placed between two outer points
    xi=ri*np.cos(ai)
    yi=ri*np.sin(ai)

    X=[]
    Y=[]
    for i in range(N):
        X+=[xo[i], xi[i]]
        Y+=[yo[i], yi[i]]

    X=X+X[:3] # duplicate first points at the end to ensure correct interpolation at the edges
    Y=Y+Y[:3] # duplicate first points at the end to ensure correct interpolation at the edges

    X=np.array(X)
    Y=np.array(Y)

    F=np.vstack((X,Y))
    s=np.arange(F.shape[1]) # parametrisation of the curve F
    F=np.transpose(F)

    Fn=cardinalCubicHermiteSpline(s,F,np.linspace(0,s.max(),10000),t=t,ends=False)
 
    fig=plt.figure()
    plt.plot(Fn[:,0],Fn[:,1],'k')
    plt.xlim((-3,3))
    plt.ylim((-3,3))
    plt.xticks([])
    plt.yticks([])            
    ax=plt.gca()
    ax.set_position([0.1,0.1,0.8,0.8])
    ax.spines['right'].set_linewidth(0)
    ax.spines['left'].set_linewidth(0)
    ax.spines['top'].set_linewidth(0)
    ax.spines['bottom'].set_linewidth(0)
    fig.set_figheight(6)
    fig.set_figwidth(6)
    plt.show()    


    
Nmax=30
RAN=np.random.randn(3*Nmax)
Randomize=widgets.ToggleButton(value=False,description='Randomize structure')
vw=widgets.IntSlider(min=10,max=Nmax,step=2,value=10,description='Vertices',continuous_update=False)
dw=widgets.FloatSlider(min=0,max=1,step=0.1,value=0.5,description='Distance',continuous_update=False,readout_format='.1f')
tw=widgets.FloatSlider(min=0.1,max=1,step=0.1,value=0.5,description='Tension',continuous_update=False,readout_format='.1f')
widget=interactive(contourgeneratorinteractive,R=Randomize,v=vw,d=dw,t=tw);





### GENERATOR




box_layout = widgets.Layout(display='flex',
                            flex_flow='row',
                            width='80%',
                           #border='solid',
                           flex_wrap='wrap')
datavertices = [str(v) for v in np.arange(10,31,2)]
datadistances = [str(round(v,1)) for v in np.arange(0,1.1,0.1)]
datatensions = [str(round(v,1)) for v in np.arange(0,1.1,0.1)]
firstlabel=widgets.Label('Vertices (v): ')
first = widgets.Box(children=[widgets.Checkbox(value=False, description=label, indent=False, layout=widgets.Layout(flex='auto', justify_content='flex-start')) for label in datavertices], layout=box_layout)
secondlabel=widgets.Label('Distance (d): ')
second = widgets.Box(children=[widgets.Checkbox(value=False, description=label, indent=False, layout=widgets.Layout(flex='auto', justify_content='flex-start')) for label in datadistances], layout=box_layout)
thirdlabel=widgets.Label('Tension (t): ')
third = widgets.Box(children=[widgets.Checkbox(value=False, description=label, indent=False, layout=widgets.Layout(flex='auto', justify_content='flex-start')) for label in datatensions], layout=box_layout)
generatebutton = widgets.Button(
    description='Generate set',
    disabled=False,
    button_style='', # 'success', 'info', 'warning', 'danger' or ''
    tooltip='Click me',
    icon='check' # (FontAwesome names without the `fa-` prefix)
)


v_all=np.array([int(v) for v in datavertices])
d_all=np.array([float(v) for v in datadistances])
t_all=np.array([float(v) for v in datatensions])

def generateset(b):
    
    k=1
    while os.path.exists('stimulus_set_output'+str(k)):
        k+=1
    outputdir='stimulus_set_output'+str(k)
    os.mkdir(outputdir)
    
    v_sel=v_all[[el.value for el in first.children]]
    d_sel=d_all[[el.value for el in second.children]]
    t_sel=t_all[[el.value for el in third.children]]

    for v in v_sel:
    
        N=int(v/2)
        RAN=np.random.randn(3*N)
        
        for d in d_sel:
            
            D=d*1.3

            # point radius tolerance, points are placed at the ring ± a radius tolerance following a normal with std=er:
            er=0.2*D 
            # point angle tolerance, points are uniformly distributed along the circle ± an angle tolerance following a normal with std=ea:
            ea=0.2*2*np.pi/float(N) 

            # outer ring
            ro=(1.5+D/2.)+er*RAN[:N] # outer radiuses
            #ao=2*np.pi*np.random.rand()+np.linspace(0,2*np.pi*(N-1)/float(N),N)+ea*np.random.randn(N)
            ao=np.linspace(0,2*np.pi*(N-1)/float(N),N)+ea*RAN[N:2*N] # outer angles
            ao.sort()
            xo=ro*np.cos(ao)
            yo=ro*np.sin(ao)

            # inner ring
            ri=(1.5-d/2)+er*RAN[2*N:3*N] # inner radiuses
            aux=np.array(list(ao[1:])+[ao[0]+2*np.pi])
            ai=(aux+ao)/2 # inner angles = each point is placed between two outer points
            xi=ri*np.cos(ai)
            yi=ri*np.sin(ai)

            X=[]
            Y=[]
            for i in range(N):
                X+=[xo[i], xi[i]]
                Y+=[yo[i], yi[i]]

            X=X+X[:3] # duplicate first points at the end to ensure correct interpolation at the edges
            Y=Y+Y[:3] # duplicate first points at the end to ensure correct interpolation at the edges

            X=np.array(X)
            Y=np.array(Y)

            F=np.vstack((X,Y))
            s=np.arange(F.shape[1]) # parametrisation of the curve F
            F=np.transpose(F)
            
            for t in t_sel:

                Fn=cardinalCubicHermiteSpline(s,F,np.linspace(0,s.max(),10000),t=t,ends=False)

                fig=plt.figure()
                plt.plot(Fn[:,0],Fn[:,1],'k')
                plt.xlim((-3,3))
                plt.ylim((-3,3))
                plt.xticks([])
                plt.yticks([])            
                ax=plt.gca()
                ax.set_position([0.1,0.1,0.8,0.8])
                ax.spines['right'].set_linewidth(0)
                ax.spines['left'].set_linewidth(0)
                ax.spines['top'].set_linewidth(0)
                ax.spines['bottom'].set_linewidth(0)
                fig.set_figheight(6)
                fig.set_figwidth(6)
                fig.savefig(outputdir+'/stimulus_v'+str(v)+'_d{:0>2.0f}_t{:0>2.0f}'.format(10*d,10*t)+'.png',dpi=200)
                plt.close()
    
    print('Done')
    
generatebutton.on_click(generateset)

generator=widgets.VBox([firstlabel,first,secondlabel,second,thirdlabel,third,generatebutton])




    
      
