# About

We aimed to use traditional computer vision methods to extract high-level features from broadcast tennis video, and use these features to analyze the player play-style pattern. We extracted features like 2D ball features based on data association method by F Yan, 3D ball trajectory using numerical optimization, player location.

## Dataset
We used matches from Australian Open 2017 videos that were publically available on the web. These are complete broadcast videos, contains gameplay, commercials, crowd, etc.

## Pre-processing Steps
Before attempting to extract the features, we need to segment out only the gameplay segments and ignore the rest, like crowd, commercials, etc. This pre-processing was done by finding the SIFT features of an initial court template and matching every frame's SIFT to get the gameplay segments. Further, stabilization of the videos is required. By computing a geometric transformation matrix from the matching SIFT features, followed by the actual transformation and interpolation, we achieve video stabilization.

## 2D ball tracking using data association method
The method starts by taking a temporal difference across frames to find ball candidates. Refinement of these ball candidates based on the size and shape of the blob is carried out. This filtering is based on thresholding the number of pixels on each blob, and fitting an ellipse in each of the blobs and filtering out high eccentricity ellipses and non-elliptical blobs. These blobs are aggregated across frames by fitting a constant acceleration model to form tracklets. Later these tracklets are optimized to grow to find better set of tracklets. Finally, these tracklets are grouped to get the ball trajectory in the video sequence. For more details see F Yans [thesis](https://pdfs.semanticscholar.org/d135/b747e99e6a06f5ecac5462b53c1b7bd259e2.pdf) and my implementaiton of it on this [page](https://researchweb.iiit.ac.in/~vishal.tiwari/ball_tracking_BTV.html).

## 3D ball Trajectory
As the videos are monocular, its impossible to find the 3D trajectory of the ball, unless we include the ball physics. We use the momentum equation to get three coupled ordinary differential equations are shown below, whose numerical solution for a set of translation velocity, angular velocity, the coefficient of restitution and coefficient of friction will give a 3D trajectory. Boundry conditions have to be taken care of, for example when the ball hits the court, the angular and translational velocities need to be changed, depending on the coefficient of friction and restitution.

![Diff Equation](https://i.imgur.com/fV4zaAT.jpg)

Now to find the 3D trajectory corresponding to a 2D path, a least squared loss function is defined between a projected initialized trajectory to 2D(using a camera calibration matrix) and the corresponding 2D path. The parameters that minimize our loss function gives us the 3D trajectory of the ball. This objective function is minimized using the interior point method in Matlab. Seen below is such a trajectory. The first pair shows the initial trajectory, which is projected to 2D (red) and the actual 2D path is shown in green, while in the 2nd row, displayed is the optimized trajectory. Note the bump at the bottom; it's because of time discretization, which defines the resolution of the trajectory.

![Init](https://i.imgur.com/5BsJDIZ.png)
![Final Trajectory](https://i.imgur.com/phZR1NP.png)

