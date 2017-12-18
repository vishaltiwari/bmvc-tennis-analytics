## About

We aimed to use traditional computer vision methods to extract high-level features from broadcast tennis video, and use these features to analyze the player play-style pattern. We extracted features like 2D ball features based on data association method by F Yan, 3D ball trajectory using numerical optimization, player location.

## Dataset
We used matches from Australian Open 2017 videos that were publically available on the web. These are complete broadcast videos, contains gameplay, commercials, crowd, etc.

## Pre-processing Steps
Before attempting to extract the features, we need to segment out only the gameplay segments and ignore the rest, like crowd, commercials, etc. This pre-processing was done by finding the SIFT features of an initial court template and matching every frame's SIFT to get the gameplay segments. Further, stabilization of the videos is required. By computing a geometric transformation matrix from the matching SIFT features, followed by the actual transformation and interpolation, we achieve video stabilization.

## 2D ball tracking using data association method
The method starts by taking a temporal difference across frames to find ball candidates. Refinement of these ball candidates based on the size and shape of the blob is carried out. This filtering is based on thresholding the number of pixels on each blob, and fitting an ellipse in each of the blobs and filtering out high eccentricity ellipses and non-elliptical blobs. These blobs are aggregated across frames by fitting a constant acceleration model to form tracklets. Later these tracklets are optimized to grow to find better set of tracklets. Finally, these tracklets are grouped to get the ball trajectory in the video sequence. For more details see F Yans [thesis](https://pdfs.semanticscholar.org/d135/b747e99e6a06f5ecac5462b53c1b7bd259e2.pdf) and my implementaiton of it on this [page](https://researchweb.iiit.ac.in/~vishal.tiwari/ball_tracking_BTV.html).
Results: We annotated a total of 1338 frames, and detecting a ball position with 5-pixel radius is marked as detection. Every frame contains a ball. Accuracy is the ratio of detections, and we get an accuracy of 0.9006. Shown below in yellow circles are ball detection from ball candidates and blacks circles are from the tracklets. 

![Ball Tracking](https://i.imgur.com/7YgmKof.png)

## 3D ball Trajectory
As the videos are monocular, its impossible to find the 3D trajectory of the ball, unless we include the ball physics. We use the momentum equation to get three coupled ordinary differential equations are shown below, whose numerical solution for a set of translation velocity, angular velocity, the coefficient of restitution and coefficient of friction will give a 3D trajectory. Boundry conditions have to be taken care of, for example when the ball hits the court, the angular and translational velocities need to be changed, depending on the coefficient of friction and restitution.

![Diff Equation](https://i.imgur.com/fV4zaAT.jpg)

Now to find the 3D trajectory corresponding to a 2D path, a least squared loss function is defined between a projected initialized trajectory to 2D(using a camera calibration matrix) and the corresponding 2D path. The parameters that minimize our loss function gives us the 3D trajectory of the ball. This objective function is minimized using the interior point method in Matlab. Seen below is such a trajectory. The first pair shows the initial trajectory, which is projected to 2D (red) and the actual 2D path is shown in green, while in the 2nd row, displayed is the optimized trajectory. Note the bump at the bottom; it's because of time discretization, which defines the resolution of the trajectory. As we didn't have actual 3D ball trajectory, the results were qualitatively evaluated based on the player's height and where he is hitting the ball.

![Init](https://i.imgur.com/5BsJDIZ.png)
![Final Trajectory](https://i.imgur.com/phZR1NP.png)

## Player Detection
One of the features of interest is the player's location in the video sequences. Our dataset is a singles tournament, so only two players, the near and the far player need to be detected. We start by extracting possible bounding boxes of players by running a pre-trained faster-RCNN for person class on each of the frames. These boxes contain players along with some of the side umpires, the match referee, etc. By using a suitable assumption that only players move during a game sequence, non-player bounding boxes can be filtered out. When the faster-RCNN fails to detect any of the players in a frame, the previously detected bounding box(with some degree of expansion) combined with foreground blobs which represent players are used to locate the location of the player in the frame. This works because players couldn't move large distances in a matter of a single frame. These foreground blobs are also used to refine the faster-RCNN results by taking the minimum bounding box of the player's foreground pixels.

The results were not satisfactory, so in the end a YOLO network was fine tuned to detect near and far-end player.(Carried out by Anurag)

![Player Tracking](https://i.imgur.com/JSZWUy5.png)

## Event Detection
Based on the y coordinates of the ball (x,y, frame) need to be classified as bounce, hit or net. This naive approach uses graphs as shown below.
Legend:
- The yellow dots show the found local maxima.
- The purple diamond shows the local minima.
- The red diamonds are the local peaks in between local maxima and local minima.
There are some issues here. As one goes from local minima to maxima (upper player hits the ball to the lower player), there are no kinks that can be found. There are no features which can be used to detect a bounce. Thus better event detection methods need to be incorporated to get better results and perform player analysis.
![Event Detection](https://i.imgur.com/WAzd0Cb.png)

## Results
From the rally in the player detection and ball tracking, we can see that the upper player is in control of the rally and making the lower player run around the court to gain a point, and he does win the point.
