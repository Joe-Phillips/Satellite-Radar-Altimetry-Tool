#----------------------------------------------------------------------
# Imports 
#----------------------------------------------------------------------

import streamlit as st
from scipy.interpolate import interp1d
import scipy.optimize as spopt
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from PIL import Image

#----------------------------------------------------------------------
# Plotting function
#----------------------------------------------------------------------

def altimetry_plot(topography,
                range_window,
                NUM_RAYS,
                ANIMATION_LENGTH,
                FPS,
                PULSE_EFFECT_DURATION,
                NOISE_PEAK,
                RAY_ANGLE_DROPOFF_MIN,
                SAT_ALTITUDE=4,
                SIMULATED_NUM_RAYS=64,
                NUM_WF_SAMPLES=128,
                ADDED_END_TIME=0.5,
                PLOT_HEIGHT=1000):

    #----------------------------------------------------------------------
    # Intake topography, normalise, and generate spline fit
    #----------------------------------------------------------------------

    topography = topography/np.max(topography)
    NUM_TOPOGRAPHY_POINTS = len(topography)
    topography_x = np.linspace(-1,1,NUM_TOPOGRAPHY_POINTS)
    spline = interp1d(topography_x, topography,fill_value="extrapolate",kind="slinear")

    #----------------------------------------------------------------------
    # Make rays of same length, originating from satellite, with uniformly spaced x
    #----------------------------------------------------------------------

    rays_x1 = np.repeat(0,SIMULATED_NUM_RAYS)
    rays_x2 = np.linspace(-1,1,SIMULATED_NUM_RAYS)
    rays_x2[np.argwhere(rays_x2==0)] = 0.00001 # adjustment so root can be found where delta x = 0
    rays_x = np.column_stack([rays_x1,rays_x2])
    rays_y1 = np.repeat(SAT_ALTITUDE,SIMULATED_NUM_RAYS)
    rays_y2 = SAT_ALTITUDE - np.sqrt(SAT_ALTITUDE**2 - np.square(rays_x2))
    rays_y = np.column_stack([rays_y1,rays_y2])

    #----------------------------------------------------------------------
    # Find bisections of rays with topography, determine ray unit vectors, and find POCA
    #----------------------------------------------------------------------

    bisections = np.full((SIMULATED_NUM_RAYS,2),np.nan)
    ray_unit_vecs = np.full((SIMULATED_NUM_RAYS,2),np.nan)
    ray_lengths = np.full(SIMULATED_NUM_RAYS,np.nan)
    prev_length = np.inf
    poca_index = 0
    for ray in range(SIMULATED_NUM_RAYS):

        # find bisections
        ray_func = interp1d(rays_x[ray], rays_y[ray],fill_value="extrapolate")
        f = lambda x: spline(x) - ray_func(x) # difference function, its zero marks the intersection
        try:
            root = spopt.root_scalar(f, x0=rays_x[ray], bracket=[-1,1]).root # find root
            bisections[ray] = [root,spline(root)]
        except:
            continue

        # determine ray unit vector
        ray_vec = np.array((rays_x[ray][0], rays_y[ray][0])) - np.array((rays_x[ray][-1], rays_y[ray][-1]))
        ray_unit_vecs[ray] = ray_vec/np.linalg.norm(ray_vec)

        # find POCA
        ray_lengths[ray] = np.linalg.norm( ray_vec - np.array(bisections[ray]) )
        if ray_lengths[ray] < prev_length:
            prev_length = ray_lengths[ray]
            poca_index = ray

    #----------------------------------------------------------------------
    # Get the reflection and return frames for each ray
    #----------------------------------------------------------------------

    # only take frames where pulses are moving, ignoring the additonal added time at the end
    ADJUSTED_ANIMATION_LENGTH = ANIMATION_LENGTH - ADDED_END_TIME
    ADJUSTED_NUM_FRAMES = int(FPS*(ADJUSTED_ANIMATION_LENGTH))

    # get reflection and return frames based on propotion of length to sat altitude
    return_frame = (ray_lengths/SAT_ALTITUDE) * ADJUSTED_NUM_FRAMES
    reflection_frame = return_frame/2
    return_frame = np.round(return_frame).astype("int")
    reflection_frame = np.round(reflection_frame).astype("int")

    #----------------------------------------------------------------------
    # Generate waveform
    #----------------------------------------------------------------------

    NUM_FRAMES = int(FPS*(ANIMATION_LENGTH))
    waveform = np.zeros(NUM_WF_SAMPLES)

    # convert return frames to samples on waveform
    FRAME_TO_SAMPLE = NUM_WF_SAMPLES/NUM_FRAMES
    return_sample = return_frame*FRAME_TO_SAMPLE

    # get the pulse drop off over time 
    PULSE_EFFECT_DURATION_SAMPLES = int(np.ceil(FPS*PULSE_EFFECT_DURATION * FRAME_TO_SAMPLE))

    if PULSE_EFFECT_DURATION_SAMPLES == 0:
        PULSE_EFFECT_DURATION_SAMPLES = 1
    pulse_drop_off_time = np.flip(np.logspace(np.log(RAY_ANGLE_DROPOFF_MIN), np.log(1), PULSE_EFFECT_DURATION_SAMPLES, base=np.exp(1)))

    # get the pulse drop off over angle 
    ODD = NUM_RAYS%2
    lhs = np.logspace(np.log(RAY_ANGLE_DROPOFF_MIN), np.log(1), SIMULATED_NUM_RAYS//2 + ODD, base=np.exp(1))
    rhs = np.flip(lhs)
    if ODD:
        pulse_drop_off_angle = np.concatenate((lhs[:-1],[1],rhs[1:]))
    else:
        pulse_drop_off_angle = np.concatenate((lhs,rhs))

    # loop through and add to waveform
    for ray in range(SIMULATED_NUM_RAYS):
        if return_sample[ray] < 0: # nan check - converted to negative overflow (I think?) in int conversion
            continue
        impact = np.zeros(NUM_WF_SAMPLES)
        pulse = pulse_drop_off_time*pulse_drop_off_angle[ray]
        try:
            impact[int(return_sample[ray]):int(return_sample[ray])+PULSE_EFFECT_DURATION_SAMPLES] = pulse
            waveform += impact
        except:
            indices_left = len(impact[int(return_sample[ray]):])
            impact[int(return_sample[ray]):] = pulse[:indices_left]
            waveform += impact

    # add noise and normalise
    noise = np.random.uniform(low=0, high=NOISE_PEAK*np.max(waveform), size=NUM_WF_SAMPLES)
    waveform += noise
    waveform = waveform/np.max(waveform)

    # convert from samples to frames for animation
    waveform = np.array([ waveform[round(frame)] for frame in np.linspace(0,NUM_WF_SAMPLES-1,NUM_FRAMES)])

    #----------------------------------------------------------------------
    # Generate location frames for animated ray pulses
    #----------------------------------------------------------------------

    # evenly sample NUM_RAY number of the simulated rays for animating
    animated_ray_indices = np.round( np.linspace(0,SIMULATED_NUM_RAYS-1,NUM_RAYS) ).astype("int")
    
    # get pulse locations
    pulse_loc = np.zeros((NUM_RAYS,NUM_FRAMES,2))

    for ray in range(NUM_RAYS):
        pulse_loc[ray,:,:] = [0,SAT_ALTITUDE]
        if reflection_frame[animated_ray_indices][ray] < 0:
            continue
        x_coords = np.linspace(0,bisections[animated_ray_indices][ray][0],reflection_frame[animated_ray_indices][ray])
        x_coords = np.concatenate((x_coords,np.flip(x_coords)))
        y_coords = np.linspace(SAT_ALTITUDE,bisections[animated_ray_indices][ray][1],reflection_frame[animated_ray_indices][ray])
        y_coords = np.concatenate((y_coords,np.flip(y_coords)))
        coords = np.dstack((x_coords,y_coords)).squeeze()

        # remove or add frames to handle rounding errors
        length_diff = return_frame[animated_ray_indices][ray] - len(coords)
        if length_diff > 0:
            coords = np.vstack((coords,np.tile([0,SAT_ALTITUDE],(length_diff,1))))
        elif length_diff < 0:
            coords = coords[:length_diff]

        pulse_loc[ray,:return_frame[animated_ray_indices][ray]] = coords

    pulse_loc = np.moveaxis(pulse_loc,1,0) # reorganise so frames are first

    #----------------------------------------------------------------------
    # Initialise figures
    #----------------------------------------------------------------------

    # get POCA on animated ray
    poca_index = animated_ray_indices[np.argmin(np.absolute(animated_ray_indices-poca_index))]

    # get ray colours, fade with angle drop off#
    ray_colours = np.array(["rgba(255, 0, 0, "+str(pulse_drop_off_angle[animated_ray_indices][ray])+")" for ray in range(NUM_RAYS)])

    # generate range window points for plot 1
    range_window_plot1 = np.array(range_window)*SAT_ALTITUDE
    range_window_points_start = np.full((NUM_RAYS,2),np.nan)
    range_window_points_end = np.full((NUM_RAYS,2),np.nan)
    for ray in range(NUM_RAYS):
        range_window_points_start[ray] = np.array((rays_x2[animated_ray_indices][ray],rays_y2[animated_ray_indices][ray])) + np.array(ray_unit_vecs[animated_ray_indices][ray]*range_window_plot1[0])
        range_window_points_end[ray] = np.array((rays_x2[animated_ray_indices][ray],rays_y2[animated_ray_indices][ray])) + np.array(ray_unit_vecs[animated_ray_indices][ray]*range_window_plot1[1])

    # generate initial elements for plot 1
    fig = make_subplots(rows=1, cols=2,column_widths=[0.65, 0.35],horizontal_spacing = 0.01)
    fig.add_trace(go.Scatter(x=topography_x,y=spline(topography_x),mode="lines",name="Topography",fill='tozeroy',showlegend=False,marker=dict(color="gray"),line_shape='spline'),row=1, col=1)
    for ray in range(NUM_RAYS): 
        fig.add_trace(go.Scatter(x=[0,bisections[animated_ray_indices][ray,0]],y=[SAT_ALTITUDE,bisections[animated_ray_indices][ray,1]],mode="lines",showlegend=False,name="Ray "+str(ray),marker=dict(color=ray_colours[ray]),line = dict(dash = 'dash')),row=1, col=1)
        fig.add_trace(go.Scatter(x=[bisections[animated_ray_indices][ray,0]],y=[bisections[animated_ray_indices][ray,1]],mode="markers",showlegend=False,name="Ray Bisect "+str(ray),marker=dict(color=ray_colours[ray],symbol="x")),row=1, col=1)
    fig.add_trace(go.Scatter(x=[bisections[poca_index,0]],y=[bisections[poca_index,1]],mode="markers",showlegend=False,name="POCA",marker=dict(color="yellow",size=12,symbol="star",line=dict(width=1.5,color='orange'))),row=1, col=1)
    fig.add_trace(go.Scatter(x=range_window_points_start[:,0],y=range_window_points_start[:,1],mode="lines",line_shape='spline',name="Range Window Start (Plot 1)",fillcolor="rgba(100, 255, 100, 0.15)",marker=dict(color="rgba(0, 255, 0, 1)"),line = dict(dash = 'dash'),showlegend=False),row=1,col=1)
    fig.add_trace(go.Scatter(x=range_window_points_end[:,0],y=range_window_points_end[:,1],mode="lines",line_shape='spline',name="Range Window End (Plot 1)",fill='tonexty',fillcolor="rgba(100, 255, 100, 0.15)",marker=dict(color="rgba(0, 255, 0, 1)"),line = dict(dash = 'dash'),showlegend=False),row=1,col=1)

    # generate range window points for plot 2
    WAVEFORM_LENGTH_WITHOUT_ADDED_TIME = NUM_FRAMES - NUM_FRAMES*(ADDED_END_TIME/ANIMATION_LENGTH)
    range_window_plot2 = np.flip(1-np.array(range_window))*WAVEFORM_LENGTH_WITHOUT_ADDED_TIME

    # generate initial elements for plot 2
    fig.add_trace(go.Scatter(x=[range_window_plot2[0],range_window_plot2[0]],y=[0,1],mode="lines",name="Range Window Start (Plot 2)",fillcolor="rgba(100, 255, 100, 0.15)",marker=dict(color="rgba(0, 255, 0, 1)"),line = dict(dash = 'dash'),showlegend=False),row=1,col=2)
    fig.add_trace(go.Scatter(x=[range_window_plot2[1],range_window_plot2[1]],y=[0,1],mode="lines",name="Range Window End (Plot 2)",fillcolor="rgba(100, 255, 100, 0.15)",fill='tonextx',marker=dict(color="rgba(0, 255, 0, 1)"),line = dict(dash = 'dash'),showlegend=False),row=1,col=2)
    fig.add_trace(go.Scatter(y=[waveform[0]],name="Waveform",fill='tozeroy',marker=dict(color="skyblue"),showlegend=False),row=1,col=2)

    # generate pulses for plot 1
    for ray in range(NUM_RAYS):
        fig.add_trace(go.Scatter(x=[0],y=[SAT_ALTITUDE],name="PULSE_"+str(ray),mode="markers",showlegend=False,marker=dict(color=ray_colours[ray])),row=1, col=1)

    #----------------------------------------------------------------------
    # Add S3 image
    #----------------------------------------------------------------------

    satImage = Image.open("s3.png")
    fig.add_layout_image(dict(source=satImage,xref="x",yref="y",x=0.03, y=0.945*SAT_ALTITUDE,sizex=0.5, sizey=0.5,xanchor="right", yanchor="bottom"))

    #----------------------------------------------------------------------
    # Make frames
    #----------------------------------------------------------------------
    
    frame_data = np.tile(go.Scatter(x=[0],y=[SAT_ALTITUDE]),(NUM_FRAMES,NUM_RAYS))
    for frame in range(NUM_FRAMES):
        for ray in range(NUM_RAYS):
            if frame >= return_frame[animated_ray_indices][ray]:
                continue
            else:
                frame_data[frame,ray] = go.Scatter(x=[pulse_loc[frame,ray,0]],y=[pulse_loc[frame,ray,1]])
 
    waveform_rw1_frame_data = np.array([go.Scatter(x=[range_window_plot2[0],range_window_plot2[0]],y=[0,1])for frame in range(NUM_FRAMES)]) # adding unchanging rw lines on waveform as otherwise they are removed on animation play
    waveform_rw2_frame_data = np.array([go.Scatter(x=[range_window_plot2[1],range_window_plot2[1]],y=[0,1])for frame in range(NUM_FRAMES)])
    waveform_frame_data = np.array([go.Scatter(y=waveform[0:frame]) for frame in range(NUM_FRAMES)])
    frame_data = np.concatenate([waveform_rw1_frame_data[:,np.newaxis],waveform_rw2_frame_data[:,np.newaxis],waveform_frame_data[:,np.newaxis], frame_data], axis=1)

    NUM_TRACES_BEFORE = 2*NUM_RAYS + 4
    frames = [dict(name = frame, data = list(frame_data[frame]), traces = np.arange(NUM_TRACES_BEFORE,NUM_RAYS+NUM_TRACES_BEFORE+3)) for frame in range(NUM_FRAMES)]

    #----------------------------------------------------------------------
    # Make sliders and buttons
    #----------------------------------------------------------------------

    updatemenus = [{
                    "buttons": [{"args": [None, {"frame": {"duration": 1000/FPS, "redraw": False},
                                                                    "fromcurrent": True, 
                                                                    "transition": {"duration":0},
                                                                    "mode":'next'}],
                                "label": "Play",
                                "method": "animate"},

                                {"args": [[None], {"frame": {"duration": 0, "redraw": False},
                                            "mode": "immediate",
                                            "transition": {"duration": 0}}],
                                "label": "Pause",
                                "method": "animate"}],

                    "direction": "left",
                    "pad": {"r": 0, "t": 50},
                    "showactive": False,
                    "type": "buttons",
                    "x": 0.5,
                    "xanchor": "right",
                    "y": 0,
                    "yanchor": "top"
                }]

    #----------------------------------------------------------------------
    # Finalise and save
    #----------------------------------------------------------------------

    fig.frames=frames
    fig.update_layout(updatemenus=updatemenus)
    fig.update_layout(xaxis2=dict(showgrid=False,visible=False,range=[-5,NUM_FRAMES]),yaxis2=dict(showgrid=False,visible=False,range=[0,1.5]))
    fig.update_layout(xaxis=dict(showgrid=False,visible=False,range=[-1,1]),yaxis=dict(showgrid=False,visible=False,range=[-0.1,4.35]))
    fig.update_layout(dict(plot_bgcolor= "rgba(0, 0, 0, 0)",paper_bgcolor= "rgba(0, 0, 0, 0)"))
    fig.update_layout(height=PLOT_HEIGHT)

    return fig

#----------------------------------------------------------------------
# Main
#----------------------------------------------------------------------

def main():

    #----------------------------------------------------------------------
    # Set up
    #----------------------------------------------------------------------

    APP_TITLE = 'Satellite Radar Altimetry Tool'

    st.set_page_config(APP_TITLE, page_icon=":satellite:", layout="wide")
    st.title(APP_TITLE)

    st.markdown(
        """
        <style>
            [data-testid=stSidebar] [data-testid=stImage]{
                text-align: center;
                display: block;
                margin-left: auto;
                margin-right: auto;
                width: 100%;
            }
        </style>
        """, unsafe_allow_html=True
    )

    st.sidebar.title(":globe_with_meridians: About")
    st.sidebar.info(
        """
        This iteractive panel web app is designed to provide an educational tool that helps visualise the complex nauances of satellite radar altimetry.
        """
    )

    st.sidebar.title(":email: Contact")
    st.sidebar.info(
        """
        Made by Joe Phillips.

        [![Repo](https://badgen.net/badge/icon/GitHub/green?icon=github&label)](https://github.com/Joe-Phillips) 
        [![Repo](https://badgen.net/badge/icon/linkedin/blue?icon=linkedin&label)](https://www.linkedin.com/in/joe-b-phillips/)
        j.phillips5@lancaster.ac.uk

        Special thanks to Dom Hardy for help with setting up Streamlit. 
        d.j.hardy@lancaster.ac.uk 
        """
    )

    st.markdown(""" <style> .font {
    font-size:20px} 
    </style> """, unsafe_allow_html=True)


    st.markdown(
        """
        ### :satellite: What is Satellite Radar Altimetry?

Satellite radar altimetry is a technique used to measure the height of surfaces from space. It works by emitting **radar pulses** :large_orange_circle: down from the satellite toward the **surface** :black_square_button: at the speed of light. These pulses bounce off the surface and return to the satellite. By measuring the time it takes for the pulses to return we can precisely calculate the distance to the surface, known as the range. Since the satellites orbit at a known altitude above the Earth, subtracting the range measurement from the satellite altitude gives the height of the surface below. By taking continuous measurements as it orbits, the satellite builds up a detailed picture of the topography of the surface. Scientists can then use this data to track changes in sea level rise, melting ice sheets, flooding, and other environmental changes.

Reflected echoes are captured in the form of a **waveform** :large_blue_square:, which records the power recieved by the altimeter over time. In general, surface elevation values extracted from the waveform are attributed to the point of closest approach (**POCA** :star:) of the surface to the satellite. These commonly correspond to a point on the foremost peak of the waveform, known as the leading edge. Although there exist many algorithms to automate this process, finding POCA and extracting associated elevation measurements from the waveform becomes more difficult over increasingly complex surfaces.

Once a pulse is emitted from the satellite, the altimeter can only measure the reflected echoes over a limited time window or **range window** :large_green_square:. If the satellite measures the echoes at the wrong time, returns can be missed. This is known as losing track. The size of this range window varies from satellite to satellite, and knowing where to place it can be a non-trivial problem, especially over complex surfaces.

        """
        )
        
    st.markdown("")
        
    st.markdown(
        """
To run the model, input a list of numbers below, representing the height of topopgraphy equidistant along the x-axis. Any list with length greater than two is allowed, and as input numbers are normalised, any number is acceptable.

First try creating an echo for a flat surface by adding **1,1** and clicking **Play** below!

        """
        )

    st.markdown(
    """
    <style>
        [data-testid=stSidebar] [data-testid=stImage]{
            text-align: center;
            display: block;
            margin-left: auto;
            margin-right: auto;
            width: 100%;
        }
    </style>
    """ , unsafe_allow_html=True
    )

    #----------------------------------------------------------------------
    # Parameters
    #----------------------------------------------------------------------

    st.markdown("""
    ### :wrench: Parameters
    
    """
    )

    topography_input = st.text_input("Topography (2 or More Comma-Seperated Numbers): ",placeholder="1,2,3,2,1")
    topography = np.array(topography_input.replace(" ","").split(","))

    try:
        topography = topography.astype("float")

        if len(topography) <= 1:
            st.write(":warning: Invalid Input Topography :warning:")
            topography = [""]

    except:
        st.write(":warning: Invalid Input Topography :warning:")

    range_window = st.slider('Range Window: ', 0.0, 1.0, (0.0, 0.4), step=0.01)
    NUM_RAYS = st.slider('Number of Rays: ', 1, 64, 16)
    NOISE_PEAK = st.slider('Noise Peak: ', 0.0, 0.1, 0.01, step=0.01)
    RAY_ANGLE_DROPOFF_MIN = st.slider('Minimum Drop-off with Ray Angle: ', 0.01, 1.0, 0.1, step=0.01)

    PULSE_EFFECT_DURATION = 0.1 #st.slider('Pulse Effect Duration (s): ', 0.0, 1.0, 0.1, step=0.01)
    ANIMATION_LENGTH = 4 #st.slider('Animation Length (s): ', 1.0, 8.0, 4.0, step=0.5)
    FPS = 30 #st.slider('FPS: ', 1, 60, 30, step=1)

    #----------------------------------------------------------------------
    # Display
    #----------------------------------------------------------------------

    st.markdown("### :chart_with_upwards_trend: Plot")

    PLOT_HEIGHT = 800
    try:
        st.plotly_chart(altimetry_plot(topography,range_window,NUM_RAYS,ANIMATION_LENGTH,FPS,PULSE_EFFECT_DURATION,NOISE_PEAK,RAY_ANGLE_DROPOFF_MIN,PLOT_HEIGHT=PLOT_HEIGHT), theme="streamlit", use_container_width=True,height=PLOT_HEIGHT)
    except:
        st.markdown("")

#----------------------------------------------------------------------
# Run
#----------------------------------------------------------------------

if __name__ == "__main__":
    main()
