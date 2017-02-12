%Store the x,y coordinates of the court lines.
row1 = [0  0 0 ; 0 1.37 0 ; 0 5.485 0 ; 0 9.6 0 ; 0 10.97 0]; % bottom baseline
row2 = [23.77  0 0; 23.77 1.37 0; 23.77 5.485 0; 23.77 9.6 0; 23.77 10.97 0]; %top baseline
row3 = [5.485 1.37 0; 5.485 5.485 0; 5.485 9.6 0]; % bottom sever line
row4 = [18.285 1.37 0; 18.285 5.485 0; 18.285 9.6 0]; % top server line
netrow = [11.885 0.456 0; 11.885 5.485 0; 11.885 10.97 0]; % net row on court
netrowheight = [11.885 0.456 0.914; 11.885 5.485 0.914; 11.885 10.97 0.914];
worldCoordinatePoints = [row1 ; row2 ; row3 ; row4 ; netrow ; netrowheight];