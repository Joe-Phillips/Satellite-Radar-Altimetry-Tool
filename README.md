# Satellite-Radar-Altimetry-Tool

This iteractive panel web app is designed to provide an educational tool that helps visualise the complex nauances of satellite radar altimetry. 

Access the Steamlit app [here](https://satellite-radar-altimetry-tool.streamlit.app/)!
:-

Made by Joe Phillips.

[![Repo](https://badgen.net/badge/icon/GitHub/green?icon=github&label)](https://github.com/Joe-Phillips) 
[![Repo](https://badgen.net/badge/icon/linkedin/blue?icon=linkedin&label)](https://www.linkedin.com/in/joe-b-phillips/)
&nbsp;✉️ j.phillips5@lancaster.ac.uk

Special thanks to Dom Hardy for help with setting up Streamlit.
d.j.hardy@lancaster.ac.uk 

### :satellite: What is Satellite Radar Altimetry?

Satellite radar altimetry is a technique used to measure the height of surfaces from space. It works by emitting **radar pulses** 🟠 down from the satellite toward the **surface** 🔲 at the speed of light. These pulses bounce off the surface and return to the satellite. By measuring the time it takes for the pulses to return we can precisely calculate the distance to the surface, known as the range. Since the satellites orbit at a known altitude above the Earth, subtracting the range measurement from the satellite altitude gives the height of the surface below. By taking continuous measurements as it orbits, the satellite builds up a detailed picture of the topography of the surface. Scientists can then use this data to track changes in sea level rise, melting ice sheets, flooding, and other environmental changes.

Reflected echoes are captured in the form of a **waveform** 🟦, which records the power recieved by the altimeter over time. In general, surface elevation values extracted from the waveform are attributed to the point of closest approach (**POCA** ⭐) of the surface to the satellite. These commonly correspond to a point on the foremost peak of the waveform, known as the leading edge. Although there exist many algorithms to automate this process, finding POCA and extracting associated elevation measurements from the waveform becomes more difficult over increasingly complex surfaces.

Once a pulse is emitted from the satellite, the altimeter can only measure the reflected echoes over a limited time window or **range window** 🟩. If the satellite measures the echoes at the wrong time, returns can be missed. This is known as losing track. The size of this range window varies from satellite to satellite, and knowing where to place it can be a non-trivial problem, especially over complex surfaces.

### :toolbox: The Model   

To help visualise the process by which satellite radar altimetry works, the model simulates a 2D case for any arbritrary, inputted topography. This is formed of a list of numbers, and represents the height of the surface, equidistant along the x-axis. 

After inputting topography, press the **PLAY** button to run the model. This emits pulses from the satellite towards the surface, decreasing exponentially in power over emission angle, which are then reflected, and subsequently captured by the satellite. On the right, the waveform records the recieved power with respect to time, with the addition of simulated Gaussian noise. For both the 2D view and the waveform, a range window is also shown, highlighted in green.
 
For given topography, any list with length greater than two is allowed, and as input numbers are normalised, any number is acceptable.
