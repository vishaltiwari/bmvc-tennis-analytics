This is a collection of 4 matlab functions that have 
been inspired by the excellent Howard Brody's book :
"Tennis science for tennis players", 
University of Pennsilvania press, 1987. (see for example http://www.cookeassociates.com/howardbrody.html)

The function "tnsstroke" simulates and plot the trajectory 
of a tennis ball given its initial position  and velocity, 
The function is able to reproduce exactly the trajectories 
shown in chapters 5 and 6 of the book.

The function "vangle" computes the vertical acceptance angle
(chapter 5) of a particular shot given the ball initial position
and the magnitude and azimut of its initial velocity.
Beware that this function works but it's currently very slow.

The function "tnsbounce" computes the linear and angular velocity
of a ball after it has bounced with a moving or static surface.
The surface could represent the court or a racket and the input 
and output trajectories could also be plotted by the function.
Most of the bouncings described in chapters 4 and 7 are fairly 
reproduced by this function.

The function "sstrat" computes the best service strategy, given
the probabilities of the strong and weak services being good and
the probabilities of winning a point if the weak and strong 
services are good (chapter 8). 

The file tennis.mdl implements the underlying simulation of a 
flying tennis ball, including spin and aerodynamics, that it's 
used by the function tnsstroke and vangle.
By changing the ball mass, radius and aerodynamics coefficients 
it can represent any flying spherical ball (i.e. a basketball).
It is not meant ot be used directly but it could after defining
appropriate initial variables in the workspace.

Each function has its help with some ready to go examples that 
you are encouraged to try. Everything runs in Matlab 5.3 as well 
as in 6.5 (and i suppose in each version between the two).

In the (not immediate) future, it may be interesting to refine 
the rotational model of the ball by including inertia moments 
and rotational aerodynamic effects, as well as a more detailed 
model of the bounce.

Giampy, Jan 2004