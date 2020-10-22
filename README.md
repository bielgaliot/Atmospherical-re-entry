# Atmospherical-re-entry
Ballistic and lifting reentry calculation of the Apollo 11 capsule, among others.

For the codes in the first reentry project folder:
  Ballistic:
  - Input is extra ballistic coefficient, aside from those given. Graphs representing different interesting values (i.e. acceleration, dynamic pressure...) as a function of altitude will be generated automatically.
  Lifting:
  -Input is extra initial flight path angle, aside from those given. Graphs representing different interesting values (i.e. acceleration, dynamic pressure...) as a function of altitude will be generated automatically.
  Mercury:
  -Requires no external input, will produce the same graphs but for the Mercury Capsule.
  
For the codes in the second reentry project folder:

They represent the attitude dynamics of the Apollo 11 capsule upon Reentry

Uncontrolled:
Outputs deceleration, velocity, angle of attack, flight path angle (...) as a function of altitude. This code does not take into account the manouver that was performed so as not to skip out of the atmosphere. Of course, the graphs show how the capsule would have shot back into space.

Controlled:
Outputs deceleration, velocity, angle of attack, flight path angle (...) as a function of altitude. This code does take into account the flipping manouver that inverted the lift of the capsule halfway through the flight, so as to safely return to earth. Hovever, two different trim points can be observed as stable in the graphs, one of them where the heatshield would not be facing windward, which would have been the death of the astronauts inside.

Zhan:
Outputs deceleration, velocity, angle of attack, flight path angle (...) as a function of altitude. This code does take into account the flipping manouver that inverted the lift of the capsule halfway through the flight, so as to safely return to earth. Moreover, it introduces the changes to the geomentry of the capsule porposed by Zhan et al. so as to change it moment coefficients and asuure the existance of a single trim point which would situate the heatshield in a favourable position.
  
  
  
