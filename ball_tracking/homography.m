function [ h wim1 ] = homography( i1, i2)


% This function estimates 2D-2D plane projective homography between two 
% perspective images using Direct Linear Transformation, RANSAC 
% and Levenberg Marquardt optimisation.

% The format for calling upon the function is as follows:
% 
% [h wim] = homography(im1, im2);
% 
% where
% 
% im1 -> 1st Image
% im2 -> 2nd Image
% h -> Returned homography matrix
% wim -> Warped version of im1 w.r.t. im2

% Code written by B S SasiKanth (bsasikanth@gmail.com) of the Indian 
% Institute of Technology (IIT), Guwahati.



% Switch off warnings for non-singular cases of RANSAC to be encountered

warning off all;

% Finding Correspondences

a = matching(i1,i2);

% Computing Homography

h = autohomest2d(a);

% Computing Left Warped Image

wim1 = (imTrans(i1,(h)));
wim1 = imresize(wim1,[size(i1,1) size(i1,2)],'bicubic');


end

%%

function [ hmatrixauto1 ] = autohomest2d( correspondence_input )
%AUTOHOMEST2D This function computes the automatic homography between two
%perspective images using RANSAC and Levenberg Marquardt optimization.

corr_req = 4;
[corr_num, tempx, tempx] = size(correspondence_input);
max_no_inliers = 0; p = 0.99; iter_no = 0;

while(1)

    iter_no = iter_no + 1;
    inlier_detail = zeros(1,corr_num);
    
    % Random number generation

    stream = RandStream('mt19937ar','seed',sum(100*clock));
    randomvector = rand(stream,1,corr_req);

    % Choosing the correspondences

    selected_correspondences = zeros(corr_req,2,2);
    selected_indices = ceil(corr_num*randomvector);

    for i=1:corr_req
    
        selected_correspondences(i,:,:) = correspondence_input(selected_indices(i),:,:);
    
    end
    
    % Computing the temporary H matrix
    
    hmatrixtemp = homest2d(selected_correspondences);
    hmatrixtempi = inv(hmatrixtemp);
    
    % Calculating the symmetric cost and checking for inliers
    
    no_inliers = 0;
    
    for i = 1:corr_num
        
        temp1 = hmatrixtemp*[correspondence_input(i,1,1); correspondence_input(i,1,2); 1];
        temp1 = temp1/temp1(3);
        
        temp2 = hmatrixtempi*[correspondence_input(i,2,1); correspondence_input(i,2,2); 1]; %#ok<MINV>
        temp2 = temp2/temp2(3);
        
        error1 = sum(abs([correspondence_input(i,2,1); correspondence_input(i,2,2); 1] - temp1));
        error2 = sum(abs([correspondence_input(i,1,1); correspondence_input(i,1,2); 1] - temp2));

        error = error1 + error2;
        
        if(error < 20)
            no_inliers = no_inliers + 1;
            inlier_detail(i) = 1;
        end
        
    end
    
    if(no_inliers > max_no_inliers)
        max_no_inliers = no_inliers;
        max_inlier_detail = inlier_detail;
        hmatrixauto = hmatrixtemp;
    end
    
    % Dynamic updatation of the no of iterations, N
    
    w = max_no_inliers/corr_num;
    N = abs((log(1-p))/(log(1-w^4)));
    
    if(N < iter_no)
        break;
    end
    
end

% Isolating inlier correspondences

inlier_corr = zeros(max_no_inliers,2,2);
k = 1;
for i=1:corr_num
    
    if(max_inlier_detail(i) == 1)
        inlier_corr(k,:,:) = correspondence_input(i,:,:);
        k = k+1;
    end
    
end

% LevMar Optimization

hvector(1:3) = hmatrixauto(1,:);
hvector(4:6) = hmatrixauto(2,:);
hvector(7:9) = hmatrixauto(3,:);

ydata = zeros(1,max_no_inliers);

hoptim = lsqcurvefit(@myfun,hvector,inlier_corr,ydata);

% Reassigning the H matrix

hmatrixauto1 = [hoptim(1), hoptim(2), hoptim(3); hoptim(4), hoptim(5), hoptim(6); hoptim(7), hoptim(8), hoptim(9)];
hmatrixauto1 = hmatrixauto1./hmatrixauto1(3,3);

% clc;
% 
% disp('The computation exited in ');
% disp(iter_no);
% disp(' iterations');

end

%%
function [ hmatrix ] = homest2d( correspondence_input )
%HOMEST2D This function computes the 2 dimensional homography for two
%different perspective images which have been inputted.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FORMAT OF THE CORRESPONDENCE MATRIX IS AS FOLLOWS:

% correspondence_input[i][j][k]

% i = correspondence_index

% j = 1 for image1
% j = 2 for image2

% k = 1 for xcoordinate
% k = 2 for ycoordinate

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[corr_num, tempx, tempx] = size(correspondence_input);

if(corr_num < 4)
    error('Invalid no. of correspondences for homography estimation');
end


% NORMALIZE THE CORRESPONDENCES 

% Find the centroid of the correspondences

xcentroid1 = 0;
ycentroid1 = 0;
xcentroid2 = 0;
ycentroid2 = 0;

for i = 1:corr_num
    
    xcentroid1 = xcentroid1 + correspondence_input(i,1,1);
    ycentroid1 = ycentroid1 + correspondence_input(i,1,2);
    xcentroid2 = xcentroid2 + correspondence_input(i,2,1);
    ycentroid2 = ycentroid2 + correspondence_input(i,2,2);
    
end

xcentroid1 = xcentroid1/corr_num;
ycentroid1 = ycentroid1/corr_num;
xcentroid2 = xcentroid2/corr_num;
ycentroid2 = ycentroid2/corr_num;

% Generating the Translation Matrices

tr1 = [1, 0, -xcentroid1; 0, 1, -ycentroid1; 0, 0, 1];
tr2 = [1, 0, -xcentroid2; 0, 1, -ycentroid2; 0, 0, 1];

% Generating Translated Correspondences

correspondence_input_tr = zeros(corr_num,2,2);

for i = 1:corr_num
    
    temp_res = tr1*[correspondence_input(i,1,1); correspondence_input(i,1,2); 1];
    correspondence_input_tr(i,1,1) = temp_res(1); 
    correspondence_input_tr(i,1,2) = temp_res(2);

    temp_res = tr2*[correspondence_input(i,2,1); correspondence_input(i,2,2); 1];
    correspondence_input_tr(i,2,1) = temp_res(1); 
    correspondence_input_tr(i,2,2) = temp_res(2);
    
end
    
% Computing the average RMS distance of translated correspondences

rms1 = 0;
rms2 = 0;

for i = 1:corr_num
    
    rms1 = rms1 + sqrt(correspondence_input_tr(i,1,1)^2 + correspondence_input_tr(i,1,2)^2);
    rms2 = rms2 + sqrt(correspondence_input_tr(i,2,1)^2 + correspondence_input_tr(i,2,2)^2);
    
end

rms1 = rms1/corr_num;
rms2 = rms2/corr_num;

% Generating the scaling matrices

scale_factor1 = sqrt(2)/rms1;
scale_factor2 = sqrt(2)/rms2;

sc1 = [scale_factor1, 0, 0; 0, scale_factor1, 0; 0, 0, 1];
sc2 = [scale_factor2, 0, 0; 0, scale_factor2, 0; 0, 0, 1];

% Normalizing correspondences to average RMS distance sqrt(2)

correspondence_input_norm = zeros(corr_num,2,2);

for i = 1:corr_num
    
    temp_res = sc1*[correspondence_input_tr(i,1,1); correspondence_input_tr(i,1,2); 1];
    correspondence_input_norm(i,1,1) = temp_res(1); 
    correspondence_input_norm(i,1,2) = temp_res(2);

    temp_res = sc2*[correspondence_input_tr(i,2,1); correspondence_input_tr(i,2,2); 1];
    correspondence_input_norm(i,2,1) = temp_res(1); 
    correspondence_input_norm(i,2,2) = temp_res(2);
    
end

% Generating Transformation matrices and required inverse for
% denormalization

T1 = sc1*tr1;
T2 = sc2*tr2;

% T2i = inv(T2);

% Generation of the A matrix for computation of the hmatrix

A = zeros(2*corr_num, 9);

for i = 1:corr_num
    
    %A(2*i-1,:) = [0, 0, 0, -correspondence_input_norm(i,1,1), -correspondence_input_norm(i,1,2), -1, correspondence_input_norm(i,2,2)*correspondence_input_norm(i,1,1), correspondence_input_norm(i,2,2)*correspondence_input_norm(i,1,2), correspondence_input_norm(i,2,2)];
    %A(2*i,:) = [correspondence_input_norm(i,1,1), correspondence_input_norm(i,1,2), 1, 0, 0, 0, -correspondence_input_norm(i,2,1)*correspondence_input_norm(i,1,1), -correspondence_input_norm(i,2,1)*correspondence_input_norm(i,1,2), -correspondence_input_norm(i,2,1)];
    A(2*i-1,:) = [correspondence_input_norm(i,1,1), correspondence_input_norm(i,1,2), 1, 0, 0, 0, -correspondence_input_norm(i,2,1)*correspondence_input_norm(i,1,1), -correspondence_input_norm(i,2,1)*correspondence_input_norm(i,1,2), -correspondence_input_norm(i,2,1)];
    A(2*i,:) = [0, 0, 0, correspondence_input_norm(i,1,1), correspondence_input_norm(i,1,2), 1, -correspondence_input_norm(i,2,2)*correspondence_input_norm(i,1,1), -correspondence_input_norm(i,2,2)*correspondence_input_norm(i,1,2), -correspondence_input_norm(i,2,2)];
    
end

% Solving the linear equation Ax = 0 using the SVD

if ((sum(sum(isnan(A))) + (sum(sum(isinf(A))) ~= 0)))
    hvector = zeros(1,9);
else
    [tempx, tempx, V] = svd(A, 0);
    hvector = V(:,9);
end

hmatrix_norm = [hvector(1), hvector(2), hvector(3); hvector(4), hvector(5), hvector(6); hvector(7), hvector(8), hvector(9)];

% Denormalizing the hmatrix_norm

hmatrix_temp = T2\hmatrix_norm;
hmatrix = hmatrix_temp*T1;


end

%%

function [ ydata ] = myfun( x, xdata )
%MYFUN LevMar cost calculator

ydata = zeros(1,length(xdata));

hmatrixtemp = [x(1), x(2), x(3); x(4), x(5), x(6); x(7), x(8), x(9)];

for i = 1:length(xdata)
        
        temp1 = hmatrixtemp*[xdata(i,1,1); xdata(i,1,2); 1];
        temp1 = temp1/temp1(3);
        
        temp2 = hmatrixtemp\[xdata(i,2,1); xdata(i,2,2); 1];
        temp2 = temp2/temp2(3);
        
        error1 = sum(abs([xdata(i,2,1); xdata(i,2,2); 1] - temp1));
        error2 = sum(abs([xdata(i,1,1); xdata(i,1,2); 1] - temp2));

        ydata(i) = error1 + error2;
        
end

end

%%

function [ a ] = matching( I1, I2 )
%MATCHING SURF matching


  % Get the Key Points
  Options.upright=true;
  Options.tresh=0.000001;
  Ipts1=OpenSurf(I1,Options);
  Ipts2=OpenSurf(I2,Options);
  
  % Put the landmark descriptors in a matrix
  D1 = reshape([Ipts1.descriptor],64,[]); 
  D2 = reshape([Ipts2.descriptor],64,[]); 

  % Find the best matches
  err=zeros(1,length(Ipts1));
  cor1=1:length(Ipts1); 
  cor2=zeros(1,length(Ipts1));
  for i=1:length(Ipts1),
      distance=sum((D2-repmat(D1(:,i),[1 length(Ipts2)])).^2,1);
      [err(i),cor2(i)]=min(distance);
  end
  
  % Sort matches on vector distance
  [tempx, ind]=sort(err); 
  cor1=cor1(ind); 
  cor2=cor2(ind);
  
  a = zeros(length(Ipts1),2,2);
  
  for i=1:length(Ipts1)
      a(i,1,1) = Ipts1(cor1(i)).x; 
      a(i,2,1) = Ipts2(cor2(i)).x;
      a(i,1,2) = Ipts1(cor1(i)).y;
      a(i,2,2) = Ipts2(cor2(i)).y;
  end

end

%%

function ipts=OpenSurf(img,Options)
% This function OPENSURF, is an implementation of SURF (Speeded Up Robust 
% Features). SURF will detect landmark points in an image, and describe
% the points by a vector which is robust against (a little bit) rotation 
% ,scaling and noise. It can be used in the same way as SIFT (Scale-invariant 
% feature transform) which is patented. Thus to align (register) two 
% or more images based on corresponding points, or make 3D reconstructions.
%
% This Matlab implementation of Surf is a direct translation of the 
% OpenSurf C# code of Chris Evans, and gives exactly the same answer. 
% Chris Evans wrote one of the best, well structured all inclusive SURF 
% implementations. On his site you can find Evaluations of OpenSURF 
% and the C# and C++ code. http://www.chrisevansdev.com/opensurf/
% Chris Evans gave me permisson to publish this code under the (Mathworks)
% BSD license.
%
% Ipts = OpenSurf(I, Options)
%
% inputs,
%   I : The 2D input image color or greyscale
%   (optional)
%   Options : A struct with options (see below)
%
% outputs,
%   Ipts : A structure with the information about all detected Landmark points
%     Ipts.x , ipts.y : The landmark position
%     Ipts.scale : The scale of the detected landmark
%     Ipts.laplacian : The laplacian of the landmark neighborhood
%     Ipts.orientation : Orientation in radians
%     Ipts.descriptor : The descriptor for corresponding point matching
%
% options,
%   Options.verbose : If set to true then useful information is 
%                     displayed (default false)
%   Options.upright : Boolean which determines if we want a non-rotation
%                       invariant result (default false)
%   Options.extended : Add extra landmark point information to the
%                   descriptor (default false)
%   Options.tresh : Hessian response treshold (default 0.0002)
%   Options.octaves : Number of octaves to analyse(default 5)
%   Options.init_sample : Initial sampling step in the image (default 2)
%   
% Example 1, Basic Surf Point Detection
% % Load image
%   I=imread('TestImages/test.png');
% % Set this option to true if you want to see more information
%   Options.verbose=false; 
% % Get the Key Points
%   Ipts=OpenSurf(I,Options);
% % Draw points on the image
%   PaintSURF(I, Ipts);
%
% Example 2, Corresponding points
% % See, example2.m
%
% Example 3, Affine registration
% % See, example3.m
%
% Function is written by D.Kroon University of Twente (July 2010)

% Add subfunctions to Matlab Search path
%functionname='OpenSurf.m';
%functiondir=which(functionname);
%functiondir=functiondir(1:end-length(functionname));
%addpath([functiondir '/SubFunctions'])
       
% Process inputs
defaultoptions=struct('tresh',0.0002,'octaves',5,'init_sample',2,'upright',false,'extended',false,'verbose',false);
if(~exist('Options','var')), 
    Options=defaultoptions; 
else
    tags = fieldnames(defaultoptions);
    for i=1:length(tags)
         if(~isfield(Options,tags{i})),  Options.(tags{i})=defaultoptions.(tags{i}); end
    end
    if(length(tags)~=length(fieldnames(Options))), 
        warning('register_volumes:unknownoption','unknown options found');
    end
end

% Create Integral Image
iimg=IntegralImage_IntegralImage(img);

% Extract the interest points
FastHessianData.thresh = Options.tresh;
FastHessianData.octaves = Options.octaves;
FastHessianData.init_sample = Options.init_sample;
FastHessianData.img = iimg;
ipts = FastHessian_getIpoints(FastHessianData,Options.verbose);

% Describe the interest points
if(~isempty(ipts))
    ipts = SurfDescriptor_DecribeInterestPoints(ipts,Options.upright, Options.extended, iimg, Options.verbose);
end

end

%%

function D=FastHessian_BuildDerivative(r,c,t,m,b)
% This function FastHessian_BuildDerivative will ..
%
% [D] = FastHessian_BuildDerivative( r,c,t,m,b )
%  
%  inputs,
%    r : 
%    c : 
%    t : 
%    m : 
%    b : 
%  
%  outputs,
%    D : 
%  
% Function is written by D.Kroon University of Twente ()

dx = (FastHessian_getResponse(m,r, c + 1, t) - FastHessian_getResponse(m,r, c - 1, t)) / 2;
dy = (FastHessian_getResponse(m,r + 1, c, t) - FastHessian_getResponse(m,r - 1, c, t)) / 2;
ds = (FastHessian_getResponse(t,r, c) - FastHessian_getResponse(b,r, c, t)) / 2;

D = [dx;dy;ds];
end

%%

function rl=FastHessian_buildResponseLayer(rl,FastHessianData)
% This function FastHessian_buildResponseLayer will ..
%
% [rl] = FastHessian_buildResponseLayer( rl,FastHessianData )
%  
%  inputs,
%    rl : 
%    FastHessianData : 
%  
%  outputs,
%    rl : 
%  
% Function is written by D.Kroon University of Twente ()

step = fix( rl.step);                      % step size for this filter
b = fix((rl.filter - 1) / 2 + 1);         % border for this filter
l = fix(rl.filter / 3);                   % lobe for this filter (filter size / 3)
w = fix(rl.filter);                       % filter size

inverse_area = 1 / double(w * w);          % normalisation factor
img=FastHessianData.img;

[ac,ar]=ndgrid(0:rl.width-1,0:rl.height-1);
ar=ar(:); ac=ac(:);

% get the image coordinates
r = int32(ar * step);
c = int32(ac * step);

% Compute response components
Dxx =   IntegralImage_BoxIntegral(r - l + 1, c - b, 2 * l - 1, w,img) - IntegralImage_BoxIntegral(r - l + 1, c - fix(l / 2), 2 * l - 1, l, img) * 3;
Dyy =   IntegralImage_BoxIntegral(r - b, c - l + 1, w, 2 * l - 1,img) - IntegralImage_BoxIntegral(r - fix(l / 2), c - l + 1, l, 2 * l - 1,img) * 3;
Dxy = + IntegralImage_BoxIntegral(r - l, c + 1, l, l,img) + IntegralImage_BoxIntegral(r + 1, c - l, l, l,img) ...
      - IntegralImage_BoxIntegral(r - l, c - l, l, l,img) - IntegralImage_BoxIntegral(r + 1, c + 1, l, l,img);

% Normalise the filter responses with respect to their size
Dxx = Dxx*inverse_area;
Dyy = Dyy*inverse_area;
Dxy = Dxy*inverse_area;

% Get the determinant of hessian response & laplacian sign
rl.responses = (Dxx .* Dyy - 0.81 * Dxy .* Dxy);
rl.laplacian = (Dxx + Dyy) >= 0;

end

%%

function responseMap=FastHessian_buildResponseMap(FastHessianData)

% Calculate responses for the first 4 FastHessianData.octaves:
% Oct1: 9,  15, 21, 27
% Oct2: 15, 27, 39, 51
% Oct3: 27, 51, 75, 99
% Oct4: 51, 99, 147,195
% Oct5: 99, 195,291,387

% Deallocate memory and clear any existing response layers
responseMap=[]; j=0;


% Get image attributes
w = (size(FastHessianData.img,2) / FastHessianData.init_sample);
h = (size(FastHessianData.img,1)/ FastHessianData.init_sample);
s = (FastHessianData.init_sample);

% Calculate approximated determinant of hessian values
if (FastHessianData.octaves >= 1)
    j=j+1; responseMap{j}=FastHessian_ResponseLayer(w,   h,   s,   9);
    j=j+1; responseMap{j}=FastHessian_ResponseLayer(w, h, s, 15);
    j=j+1; responseMap{j}=FastHessian_ResponseLayer(w, h, s, 21);
    j=j+1; responseMap{j}=FastHessian_ResponseLayer(w, h, s, 27);
end

if (FastHessianData.octaves >= 2)
    j=j+1; responseMap{j}=FastHessian_ResponseLayer(w / 2, h / 2, s * 2, 39);
    j=j+1; responseMap{j}=FastHessian_ResponseLayer(w / 2, h / 2, s * 2, 51);
end

if (FastHessianData.octaves >= 3)
    j=j+1; responseMap{j}=FastHessian_ResponseLayer(w / 4, h / 4, s * 4, 75);
    j=j+1; responseMap{j}=FastHessian_ResponseLayer(w / 4, h / 4, s * 4, 99);
end

if (FastHessianData.octaves >= 4)
    j=j+1; responseMap{j}=FastHessian_ResponseLayer(w / 8, h / 8, s * 8, 147);
    j=j+1; responseMap{j}=FastHessian_ResponseLayer(w / 8, h / 8, s * 8, 195);
end

if (FastHessianData.octaves >= 5)
    j=j+1; responseMap{j}=FastHessian_ResponseLayer(w / 16, h / 16, s * 16, 291);
    j=j+1; responseMap{j}=FastHessian_ResponseLayer(w / 16, h / 16, s * 16, 387);
end

% Extract responses from the image
for i=1:length(responseMap);
    responseMap{i}=FastHessian_buildResponseLayer(responseMap{i},FastHessianData);
end


end

%%

function ipts=FastHessian_getIpoints(FastHessianData,verbose)
% filter index map

filter_map = [0,1,2,3;
    1,3,4,5;
    3,5,6,7;
    5,7,8,9;
    7,9,10,11]+1;

np=0; ipts=struct;

% Build the response map
responseMap=FastHessian_buildResponseMap(FastHessianData);

% Find the maxima acrros scale and space
for o = 1:FastHessianData.octaves
    for i = 1:2
        b = responseMap{filter_map(o,i)};
        m = responseMap{filter_map(o,i+1)};
        t = responseMap{filter_map(o,i+2)};
        
        % loop over middle response layer at density of the most
        % sparse layer (always top), to find maxima across scale and space
        [c,r]=ndgrid(0:t.width-1,0:t.height-1);
        r=r(:); c=c(:);
        
        p=find(FastHessian_isExtremum(r, c, t, m, b,FastHessianData));
        for j=1:length(p);
            ind=p(j);
            [ipts,np]=FastHessian_interpolateExtremum(r(ind), c(ind), t, m, b, ipts,np);
        end
    end
end

% Show laplacian and response maps with found interest-points
if(verbose)
    % Show the response map
    if(verbose)
        fig_h=ceil(length(responseMap)/3);
        h=figure;  set(h,'name','Laplacian');
        for i=1:length(responseMap), 
            pic=reshape(responseMap{i}.laplacian,[responseMap{i}.width responseMap{i}.height]);
            subplot(3,fig_h,i); imshow(pic,[]); hold on;
        end
        h=figure; set(h,'name','Responses');
        h_res=zeros(1,length(responseMap));
        for i=1:length(responseMap), 
            pic=reshape(responseMap{i}.responses,[responseMap{i}.width responseMap{i}.height]);
            h_res(i)=subplot(3,fig_h,i); imshow(pic,[]); hold on;
        end
    end
    
    % Show the maximum points
    disp(['Number of interest points found ' num2str(np)]);
    scales=zeros(1,length(responseMap));
    scaley=zeros(1,length(responseMap));
    scalex=zeros(1,length(responseMap));
    for i=1:length(responseMap)
        scales(i)=responseMap{i}.filter*(2/15);
        scalex(i)=responseMap{i}.width/size(FastHessianData.img,2);
        scaley(i)=responseMap{i}.height/size(FastHessianData.img,1);
    end
    for i=1:np
        [t,ind]=min((scales-ipts(i).scale).^2);
        plot(h_res(ind),ipts(i).y*scaley(ind)+1,ipts(i).x*scalex(ind)+1,'o','color',rand(1,3));
    end
end
end

%%

function an=FastHessian_getLaplacian(a,row, column,b)
% This function FastHessian_getLaplacian will ..
%
% [an] = FastHessian_getLaplacian( a,row,column,b )
%  
%  inputs,
%    a : 
%    row : 
%    column : 
%    b : 
%  
%  outputs,
%    an : 
%  
% Function is written by D.Kroon University of Twente ()
if(nargin<4)
    scale=1;
else
    scale=fix(a.width/b.width);
end

an=a.laplacian(fix(scale*row) * a.width + fix(scale*column)+1);

end

%%

function an=FastHessian_getResponse(a,row, column,b)
% This function FastHessian_getResponse will ..
%
% [an] = FastHessian_getResponse( a,row,column,b )
%  
%  inputs,
%    a : 
%    row : 
%    column : 
%    b : 
%  
%  outputs,
%    an : 
%  
% Function is written by D.Kroon University of Twente ()
if(nargin<4)
    scale=1;
else
    scale=fix(a.width/b.width);
end

an=a.responses(fix(scale*row) * a.width + fix(scale*column)+1);
end

%%

function [ipts, np]=FastHessian_interpolateExtremum(r, c,  t,  m,  b,  ipts, np)
% This function FastHessian_interpolateExtremum will ..
%
% [ipts,np] = FastHessian_interpolateExtremum( r,c,t,m,b,ipts,np )
%  
%  inputs,
%    r : 
%    c : 
%    t : 
%    m : 
%    b : 
%    ipts : 
%    np : 
%  
%  outputs,
%    ipts : 
%    np : 
%  
% Function is written by D.Kroon University of Twente (July 2010)
D = FastHessian_BuildDerivative(r, c, t, m, b);
H = FastHessian_BuildHessian(r, c, t, m, b);

%get the offsets from the interpolation
Of = - H\D;
O=[ Of(1, 1), Of(2, 1), Of(3, 1) ];

%get the step distance between filters
filterStep = fix((m.filter - b.filter));

%If point is sufficiently close to the actual extremum
if (abs(O(1)) < 0.5 && abs(O(2)) < 0.5 && abs(O(3)) < 0.5)
    np=np+1;
    ipts(np).x = double(((c + O(1))) * t.step);
    ipts(np).y = double(((r + O(2))) * t.step);
    ipts(np).scale = double(((2/15) * (m.filter + O(3) * filterStep)));
    ipts(np).laplacian = fix(FastHessian_getLaplacian(m,r,c,t));
end
  
function D=FastHessian_BuildDerivative(r,c,t,m,b)
dx = (FastHessian_getResponse(m,r, c + 1, t) - FastHessian_getResponse(m,r, c - 1, t)) / 2;
dy = (FastHessian_getResponse(m,r + 1, c, t) - FastHessian_getResponse(m,r - 1, c, t)) / 2;
ds = (FastHessian_getResponse(t,r, c) - FastHessian_getResponse(b,r, c, t)) / 2;
D = [dx;dy;ds];
end

function H=FastHessian_BuildHessian(r, c, t, m, b)
v = FastHessian_getResponse(m, r, c, t);
dxx = FastHessian_getResponse(m,r, c + 1, t) + FastHessian_getResponse(m,r, c - 1, t) - 2 * v;
dyy = FastHessian_getResponse(m,r + 1, c, t) + FastHessian_getResponse(m,r - 1, c, t) - 2 * v;
dss = FastHessian_getResponse(t,r, c) + FastHessian_getResponse(b,r, c, t) - 2 * v;
dxy = (FastHessian_getResponse(m,r + 1, c + 1, t) - FastHessian_getResponse(m,r + 1, c - 1, t) - FastHessian_getResponse(m,r - 1, c + 1, t) + FastHessian_getResponse(m,r - 1, c - 1, t)) / 4;
dxs = (FastHessian_getResponse(t,r, c + 1) - FastHessian_getResponse(t,r, c - 1) - FastHessian_getResponse(b,r, c + 1, t) + FastHessian_getResponse(b,r, c - 1, t)) / 4;
dys = (FastHessian_getResponse(t,r + 1, c) - FastHessian_getResponse(t,r - 1, c) - FastHessian_getResponse(b,r + 1, c, t) + FastHessian_getResponse(b,r - 1, c, t)) / 4;

H = zeros(3,3);
H(1, 1) = dxx;
H(1, 2) = dxy;
H(1, 3) = dxs;
H(2, 1) = dxy;
H(2, 2) = dyy;
H(2, 3) = dys;
H(3, 1) = dxs;
H(3, 2) = dys;
H(3, 3) = dss;

end
end

%%

function an=FastHessian_isExtremum(r, c,  t,  m,  b,FastHessianData)
% This function FastHessian_isExtremum will ..
%
% [an] = FastHessian_isExtremum( r,c,t,m,b,FastHessianData )
%  
%  inputs,
%    r : 
%    c : 
%    t : 
%    m : 
%    b : 
%    FastHessianData : 
%  
%  outputs,
%    an : 
%  
% Function is written by D.Kroon University of Twente (July 2010)

% bounds check
layerBorder = fix((t.filter + 1) / (2 * t.step));
bound_check_fail=(r <= layerBorder | r >= t.height - layerBorder | c <= layerBorder | c >= t.width - layerBorder);

% check the candidate point in the middle layer is above thresh 
candidate = FastHessian_getResponse(m,r,c,t);
treshold_fail=candidate < FastHessianData.thresh;

an=(~bound_check_fail)&(~treshold_fail);
for rr = -1:1
    for  cc = -1:1
          %  if any response in 3x3x3 is greater then the candidate is not a maximum
          check1=FastHessian_getResponse(t,r + rr, c + cc, t) >= candidate;
          check2=FastHessian_getResponse(m,r + rr, c + cc, t) >= candidate;
          check3=FastHessian_getResponse(b,r + rr, c + cc, t) >= candidate;
          check4=(rr ~= 0 || cc ~= 0);
          an3 = ~(check1 | (check4 & check2) | check3);
          an=an&an3;
    end
end

function an=FastHessian_getResponse(a,row, column,b)
scale=fix(a.width/b.width);
% Clamp to boundary 
% (The orignal C# code, doesn't contain this boundary clamp because if you 
% process one coordinate at the time you already returned on the boundary check)
index=fix(scale*row) * a.width + fix(scale*column)+1;
index(index<1)=1; index(index>length(a.responses))=length(a.responses);
an=a.responses(index);
end
end

%%

function ResponseLayerData=FastHessian_ResponseLayer(width, height, step, filter)
% This function FastHessian_ResponseLayer will ..
%
% [ResponseLayerData] = FastHessian_ResponseLayer( width,height,step,filter )
%  
%  inputs,
%    width : 
%    height : 
%    step : 
%    filter : 
%  
%  outputs,
%    ResponseLayerData : 
%  
% Function is written by D.Kroon University of Twente (July 2010)
width=floor(width);
height=floor(height);
step=floor(step);
filter=floor(filter);

ResponseLayerData.width = width;
ResponseLayerData.height = height;
ResponseLayerData.step = step;
ResponseLayerData.filter = filter;

ResponseLayerData.responses = zeros(width * height,1);
ResponseLayerData.laplacian = zeros(width * height,1);
end

%%

function an=IntegralImage_BoxIntegral(row, col, rows,cols,img)

% Get integer coordinates
row=fix(row);
col=fix(col);
rows=fix(rows);
cols=fix(cols);

% Get the corner coordinates of the box integral
r1 = min(row, size(img,1));
c1 = min(col, size(img,2));
r2 = min(row + rows, size(img,1));
c2 = min(col + cols, size(img,2));

% Get the values at the cornes of the box integral (fast 1D index look up)
sx=size(img,1);
A = img(max(r1+(c1-1)*sx,1));
B = img(max(r1+(c2-1)*sx,1));
C = img(max(r2+(c1-1)*sx,1));
D = img(max(r2+(c2-1)*sx,1));

% If coordinates are outside at the top or left, the value must be zero
A((r1<1)|(c1<1))=0;
B((r1<1)|(c2<1))=0;
C((r2<1)|(c1<1))=0;
D((r2<1)|(c2<1))=0;

% Minimum value of the integral is zero
an=max(0, A - B - C + D);


end

%%

function an=IntegralImage_HaarX(row, column, size, img)
% This function IntegralImage_HaarX will ..
%
% [an] = IntegralImage_HaarX( row,column,size,img )
%  
%  inputs,
%    row : 
%    column : 
%    size : 
%    img : 
%  
%  outputs,
%    an : 
%  
% Function is written by D.Kroon University of Twente (July 2010)
an= IntegralImage_BoxIntegral(row - size / 2, column, size, size / 2, img) - IntegralImage_BoxIntegral(row - size / 2, column - size / 2, size, size / 2, img);
end

%%

function an=IntegralImage_HaarY(row, column, size, img)
% This function IntegralImage_HaarY will ..
%
% [an] = IntegralImage_HaarY( row,column,size,img )
%  
%  inputs,
%    row : 
%    column : 
%    size : 
%    img : The integral image
%  
%  outputs,
%    an : The haar response in y-direction
%  
% Function is written by D.Kroon University of Twente (July 2010)
an= IntegralImage_BoxIntegral(row, column - size / 2, size / 2, size , img) - IntegralImage_BoxIntegral(row - size / 2, column - size / 2, size / 2, size , img);

end

%%

function pic=IntegralImage_IntegralImage(I)
% This function IntegralImage_IntegralImage will ..
%
% J = IntegralImage_IntegralImage( I )
%  
%  inputs,
%    I : An 2D image color or greyscale
%  
%  outputs,
%    J : The integral image 
%  
% Function is written by D.Kroon University of Twente (July 2010)

% Convert Image to double
switch(class(I));
    case 'uint8'
        I=double(I)/255;
    case 'uint16'
        I=double(I)/65535;
    case 'int8'
        I=(double(I)+128)/255;
    case 'int16'
        I=(double(I)+32768)/65535;
    otherwise
        I=double(I);
end

% Convert Image to greyscale
if(size(I,3)==3)
	cR = .2989; cG = .5870; cB = .1140;
	I=I(:,:,1)*cR+I(:,:,2)*cG+I(:,:,3)*cB;
end

% Make the integral image
pic = cumsum(cumsum(I,1),2);
end

%%

function PaintSURF(I, ipts)
% This function PaintSURF will display the image with the found  Interest points
%
% [] = PaintSURF( img,ipts )
%  
%  inputs,
%    img : Image 2D color or greyscale
%    ipts : The interest points
%  
% Function is written by D.Kroon University of Twente (July 2010)

% Convert Image to double
switch(class(I));
    case 'uint8'
        I=double(I)/255;
    case 'uint16'
        I=double(I)/65535;
    case 'int8'
        I=(double(I)+128)/255;
    case 'int16'
        I=(double(I)+32768)/65535;
    otherwise
        I=double(I);
end

figure, imshow(I), hold on;
if (isempty(fields(ipts))), return; end
for i=1:length(ipts)
   ip=ipts(i);
   
   S = 2 * fix(2.5 * ip.scale);
   R = fix(S / 2);

   pt =  [(ip.x), (ip.y)];
   ptR = [(R * cos(ip.orientation)), (R * sin(ip.orientation))];

   if(ip.laplacian >0), myPen =[0 0 1]; else myPen =[1 0 0]; end
   
   rectangle('Curvature', [1 1],'Position', [pt(1)-R, pt(2)-R, S, S],'EdgeColor',myPen);
   
    plot([pt(1), pt(1)+ptR(1)]+1,[pt(2), pt(2)+ptR(2)]+1,'g');
end
end

%%

function ipts = SurfDescriptor_DecribeInterestPoints(ipts, upright, extended, img, verbose)
% This function SurfDescriptor_DecribeInterestPoints will ..
%
% [ipts] = SurfDescriptor_DecribeInterestPoints( ipts,upright,extended,img )
%  
%  inputs,
%    ipts : Interest Points (x,y,scale)
%    bUpright : If true not rotation invariant descriptor
%    bExtended :  If true make a 128 values descriptor
%    img : Integral image
%    verbose : If true show useful information
%  
%  outputs,
%    ipts :  Interest Points (x,y,orientation,descriptor)
%  
% Function is written by D.Kroon University of Twente (July 2010)
if (isempty(fields(ipts))), return; end

if(verbose), h_ang=figure; drawnow, set(h_ang,'name','Angles'); else h_ang=[]; end
if(verbose), h_des=figure; drawnow, set(h_des,'name','Aligned Descriptor XY'); end
   
for i=1:length(ipts)
   % Display only information about the first 40 points
   if(i>40), verbose=false; end
   
   ip=ipts(i);
   % determine descriptor size
   if (extended), ip.descriptorLength = 128; else ip.descriptorLength = 64; end

   % Get the orientation
   if(verbose), figure(h_ang), subplot(5,8,i), end
   ip.orientation=SurfDescriptor_GetOrientation(ip,img,verbose);

   % Extract SURF descriptor
   if(verbose), figure(h_des), subplot(10,4,i), end
   ip.descriptor=SurfDescriptor_GetDescriptor(ip, upright, extended, img,verbose);
   
   ipts(i).orientation=ip.orientation;
   ipts(i).descriptor=ip.descriptor;
end

if(~isempty(h_ang)), figure(h_ang), colormap(jet); end
end

%%

function descriptor=SurfDescriptor_GetDescriptor(ip, bUpright, bExtended, img, verbose)
% This function SurfDescriptor_GetDescriptor will ..
%
% [descriptor] = SurfDescriptor_GetDescriptor( ip,bUpright,bExtended,img )
%  
%  inputs,
%    ip : Interest Point (x,y,scale, orientation)
%    bUpright : If true not rotation invariant descriptor
%    bExtended :  If true make a 128 values descriptor
%    img : Integral image
%    verbose : If true show additional information
%  
%  outputs,
%    descriptor : Descriptor of interest point length 64 or 128 (extended)  
%  
% Function is written by D.Kroon University of Twente (July 2010)

% Get rounded InterestPoint data
X = round(ip.x);
Y = round(ip.y);
S = round(ip.scale);

if (bUpright)
    co = 1;
    si = 0;
else
    co = cos(ip.orientation);
    si = sin(ip.orientation);
end

% Basis coordinates of samples, if coordinate 0,0, and scale 1
[lb,kb]=ndgrid(-4:4,-4:4); lb=lb(:); kb=kb(:);

%Calculate descriptor for this interest point
[jl,il]=ndgrid(0:3,0:3); il=il(:)'; jl=jl(:)';

ix = (il*5-8);
jx = (jl*5-8);

% 2D matrices instead of double for-loops, il, jl
cx=length(lb); cy=length(ix);
lb=repmat(lb,[1 cy]); lb=lb(:);
kb=repmat(kb,[1 cy]); kb=kb(:);
ix=repmat(ix,[cx 1]); ix=ix(:);
jx=repmat(jx,[cx 1]); jx=jx(:);

% Coordinates of samples (not rotated)
l=lb+jx; k=kb+ix;

%Get coords of sample point on the rotated axis
sample_x = round(X + (-l * S * si + k * S * co)); 
sample_y = round(Y + (l * S * co + k * S * si));
                
%Get the gaussian weighted x and y responses
xs = round(X + (-(jx+1) * S * si + (ix+1) * S * co));
ys = round(Y + ((jx+1) * S * co + (ix+1) * S * si));

gauss_s1 = SurfDescriptor_Gaussian(xs - sample_x, ys - sample_y, 2.5 * S);
rx = IntegralImage_HaarX(sample_y, sample_x, 2 * S,img);
ry = IntegralImage_HaarY(sample_y, sample_x, 2 * S,img);
                
%Get the gaussian weighted x and y responses on the aligned axis
rrx = gauss_s1 .* (-rx * si + ry * co);  rrx=reshape(rrx,cx,cy);
rry = gauss_s1 .* ( rx * co + ry * si);  rry=reshape(rry,cx,cy);
        
% Get the gaussian scaling
cx = -0.5 + il + 1; cy = -0.5 + jl + 1;
gauss_s2 = SurfDescriptor_Gaussian(cx - 2, cy - 2, 1.5);

if (bExtended)
    % split x responses for different signs of y
    check=rry >= 0; rrx_p=rrx.*check;  rrx_n=rrx.*(~check);
    
    dx = sum(rrx_p); mdx = sum(abs(rrx_p),1);
    dx_yn = sum(rrx_n); mdx_yn = sum(abs(rrx_n),1);

    % split y responses for different signs of x
    check=(rrx >= 0); rry_p=rry.*check; rry_n=rry.*(~check);
    dy = sum(rry_p,1); 
    mdy = sum(abs(rry_p),1);
    dy_xn = sum(rry_n,1); 
    mdy_xn =  sum(abs(rry_n),1);
else
    dx = sum(rrx,1);
    dy = sum(rry,1);
    mdx = sum(abs(rrx),1);
    mdy = sum(abs(rry),1);
    dx_yn = 0; mdx_yn = 0; 
    dy_xn = 0; mdy_xn = 0;
end

if (bExtended)
    descriptor=[dx;dy;mdx;mdy;dx_yn;dy_xn;mdx_yn;mdy_xn].* repmat(gauss_s2,[8 1]);
else
    descriptor=[dx;dy;mdx;mdy].* repmat(gauss_s2,[4 1]);
end
  
len = sum((dx.^2 + dy.^2 + mdx.^2 + mdy.^2 + dx_yn + dy_xn + mdx_yn + mdy_xn) .* gauss_s2.^2);

%Convert to Unit Vector
descriptor= descriptor(:) / sqrt(len);
if(verbose)
    for i=1:size(rrx,2)
        p1=reshape(rrx(:,i),[9,9]);
        p2=reshape(rry(:,i),[9,9]);
        p=[p1;ones(1,9)*0.02;p2];
        if(i==1)
            pic=p;
        else
            pic=[pic ones(19,1)*0.02 p];
        end
    end
    imshow(pic,[]);
end

function an= SurfDescriptor_Gaussian(x, y, sig)
an = 1 / (2 * pi * sig^2) .* exp(-(x.^2 + y.^2) / (2 * sig^2));

end
end

%%

function orientation=SurfDescriptor_GetOrientation(ip,img,verbose)
% This function SurfDescriptor_GetOrientation will ..
%
% [orientation] = SurfDescriptor_GetOrientation( ip,img )
%  
%  inputs,
%    ip :  InterestPoint data, (x,y,scale)
%    img : Integral Image
%    verbose : If true, show additional information
%  
%  outputs,
%    orientation : Orientation of intereset point (radians)
%  
% Function is written by D.Kroon University of Twente (July 2010)

gauss25 = [0.02350693969273 0.01849121369071 0.01239503121241 0.00708015417522 0.00344628101733 0.00142945847484 0.00050524879060;
           0.02169964028389 0.01706954162243 0.01144205592615 0.00653580605408 0.00318131834134 0.00131955648461 0.00046640341759;
           0.01706954162243 0.01342737701584 0.00900063997939 0.00514124713667 0.00250251364222 0.00103799989504 0.00036688592278;
           0.01144205592615 0.00900063997939 0.00603330940534 0.00344628101733 0.00167748505986 0.00069579213743 0.00024593098864;
           0.00653580605408 0.00514124713667 0.00344628101733 0.00196854695367 0.00095819467066 0.00039744277546 0.00014047800980;
           0.00318131834134 0.00250251364222 0.00167748505986 0.00095819467066 0.00046640341759 0.00019345616757 0.00006837798818;
           0.00131955648461 0.00103799989504 0.00069579213743 0.00039744277546 0.00019345616757 0.00008024231247 0.00002836202103];
gauss25=gauss25(:);


% Get rounded InterestPoint data
X = round(ip.x);
Y = round(ip.y);
S = round(ip.scale);

% calculate haar responses for points within radius of 6*scale
[j,i]=ndgrid(-6:6,-6:6);
j=j(:); i=i(:); check=(i.^2 + j.^2 < 36); j=j(check); i=i(check);

% Get gaussian filter (by mirroring gauss25)
id = [ 6, 5, 4, 3, 2, 1, 0, 1, 2, 3, 4, 5, 6 ];
gauss = gauss25(id(i + 6 + 1) + id(j + 6 + 1) *7+1);

resX = gauss .* IntegralImage_HaarX(Y + j * S, X + i * S, 4 * S, img);
resY = gauss .* IntegralImage_HaarY(Y + j * S, X + i * S, 4 * S, img);
Ang =  mod(atan2(resY, resX),2*pi);

% loop slides pi/3 window around feature point
ang1 = 0:0.15:(2 * pi);
ang2 = mod(ang1+pi/3,2*pi);

% Repmat is used to check for all angles (x direction) and 
% all responses (y direction) without for-loops.
cx=length(Ang); cy=length(ang1);
ang1=repmat(ang1,[cx 1]);
ang2=repmat(ang2,[cx 1]);
Ang =repmat(Ang,[1 cy]);
resX =repmat(resX,[1 cy]);
resY =repmat(resY,[1 cy]);

% determine whether the point is within the window
check1= (ang1 < ang2) & (ang1 < Ang) & (Ang < ang2);
check2= (ang2 < ang1) & ( ((Ang > 0) & (Ang < ang2)) | ((Ang > ang1) & (Ang < pi)) );
check=check1|check2;

sumX =  sum(resX.*check,1);
sumY =  sum(resY.*check,1);

% Find the most dominant direction
R=sumX.^2+ sumY.^2;
[t,ind]=max(R);
orientation =  mod(atan2(sumY(ind), sumX(ind)),2*pi);

if(verbose)
    pica=zeros(13,13); 
    pica(i+7+(j+6)*13)=Ang(:,1);
    imshow(pica,[0 2*pi]); 
end


end


%%

% IMTRANS - Homogeneous transformation of an image.
%
% Applies a geometric transform to an image
%
%  [newim, newT] = imTrans(im, T, region, sze);
%
%  Arguments: 
%        im     - The image to be transformed.
%        T      - The 3x3 homogeneous transformation matrix.
%        region - An optional 4 element vector specifying 
%                 [minrow maxrow mincol maxcol] to transform.
%                 This defaults to the whole image if you omit it
%                 or specify it as an empty array [].
%        sze    - An optional desired size of the transformed image
%                 (this is the maximum No of rows or columns).
%                 This defaults to the maximum of the rows and columns
%                 of the original image.
%
%  Returns:
%        newim  - The transformed image.
%        newT   - The transformation matrix that relates transformed image
%                 coordinates to the reference coordinates for use in a
%                 function such as DIGIPLANE.
%
%  The region argument is used when one is inverting a perspective
%  transformation of a plane and the vanishing line of the plane lies within
%  the image.  Attempts to transform any part of the vanishing line will
%  position you at infinity.  Accordingly one should specify a region that
%  excludes any part of the vanishing line.
%
%  The sze parameter is optionally used to control the size of the
%  output image.  When inverting a perpective or affine transformation
%  the scale parameter is unknown/arbitrary, and without specifying
%  it explicitly the transformed image can end up being very small 
%  or very large.
%
%  Problems: If your transformed image ends up as being two small bits of
%  image separated by a large black area then the chances are that you have
%  included the vanishing line of the plane within the specified region to
%  transform.  If your image degenerates to a very thin triangular shape
%  part of your region is probably very close to the vanishing line of the
%  plane.

% Copyright (c) 2000-2005 Peter Kovesi
% School of Computer Science & Software Engineering
% The University of Western Australia
% http://www.csse.uwa.edu.au/
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in 
% all copies or substantial portions of the Software.
%
% The Software is provided "as is", without warranty of any kind.

% April 2000 - original version.
% July 2001  - transformation of region boundaries corrected.

function [newim, newT] = imTrans(im, T, region, sze);

if isa(im,'uint8')
    im = double(im);  % Make sure image is double     
end

% Set up default region and transformed image size values
if ndims(im) == 3
    [rows cols depth] = size(im);
else
    [rows cols] = size(im);
    depth = 1;
end

if nargin == 2
    region = [1 rows 1 cols];
    sze = max([rows cols]);
elseif nargin == 3    
    sze = max([rows cols]);
end

if isempty(region)
    region = [1 rows 1 cols];
end

	
threeD = (ndims(im)==3);  % A colour image
if threeD    % Transform red, green, blue components separately
    im = im/255;  
    [r, newT] = transformImage(im(:,:,1), T, region, sze);
    [g, newT] = transformImage(im(:,:,2), T, region, sze);
    [b, newT] = transformImage(im(:,:,3), T, region, sze);
    
    newim = repmat(uint8(0),[size(r),3]);
    newim(:,:,1) = uint8(round(r*255));
    newim(:,:,2) = uint8(round(g*255));
    newim(:,:,3) = uint8(round(b*255));
    
else                % Assume the image is greyscale
    [newim, newT] = transformImage(im, T, region, sze);
end
end

%------------------------------------------------------------

% The internal function that does all the work

function [newim, newT] = transformImage(im, T, region, sze);

[rows, cols] = size(im);

if 0
% Determine default parameters if needed
if nargin == 2
  region = [1 rows 1 cols];
  sze = max(rows,cols);
elseif nargin == 3
  sze = max(rows,cols);
elseif nargin ~= 4
  error('Incorrect arguments to imtrans');
end
end
% Cut the image down to the specified region
%if nargin == 3 | nargin == 4
    im = im(region(1):region(2), region(3):region(4));
    [rows, cols] = size(im);
%end

% Find where corners go - this sets the bounds on the final image
B = bounds(T,region);
nrows = B(2) - B(1);
ncols = B(4) - B(3);

% Determine any rescaling needed
s = sze/max(nrows,ncols);

S = [s 0 0        % Scaling matrix
     0 s 0
     0 0 1];

T = S*T;
Tinv = inv(T);

% Recalculate the bounds of the new (scaled) image to be generated
B = bounds(T,region);
nrows = B(2) - B(1);
ncols = B(4) - B(3);

% Construct a transformation matrix that relates transformed image
% coordinates to the reference coordinates for use in a function such as
% DIGIPLANE.  This transformation is just an inverse of a scaling and
% origin shift. 
newT=inv(S - [0 0 B(3); 0 0 B(1); 0 0 0]);

% Set things up for the image transformation.
newim = zeros(nrows,ncols);
[xi,yi] = meshgrid(1:ncols,1:nrows);    % All possible xy coords in the image.

% Transform these xy coords to determine where to interpolate values
% from. Note we have to work relative to x=B(3) and y=B(1).
sxy = homoTrans(Tinv, [xi(:)'+B(3) ; yi(:)'+B(1) ; ones(1,ncols*nrows)]);
xi = reshape(sxy(1,:),nrows,ncols);
yi = reshape(sxy(2,:),nrows,ncols);

[x,y] = meshgrid(1:cols,1:rows);
x = x+region(3)-1; % Offset x and y relative to region origin.
y = y+region(1)-1; 
newim = interp2(x,y,double(im),xi,yi); % Interpolate values from source image.


% Plot bounding region
%P = [region(3) region(4) region(4) region(3)
%     region(1) region(1) region(2) region(2)
%      1    1    1    1   ];
%B = round(homoTrans(T,P));
%Bx = B(1,:);
%By = B(2,:);
%Bx = Bx-min(Bx); Bx(5)=Bx(1);
%By = By-min(By); By(5)=By(1);
%show(newim,2), axis xy
%line(Bx,By,'Color',[1 0 0],'LineWidth',2);
% end plot bounding region

end


%---------------------------------------------------------------------
%
% Internal function to find where the corners of a region, R
% defined by [minrow maxrow mincol maxcol] are transformed to 
% by transform T and returns the bounds, B in the form 
% [minrow maxrow mincol maxcol]

function B = bounds(T, R)

P = [R(3) R(4) R(4) R(3)      % homogeneous coords of region corners
     R(1) R(1) R(2) R(2)
      1    1    1    1   ];
     
PT = round(homoTrans(T,P)); 

B = [min(PT(2,:)) max(PT(2,:)) min(PT(1,:)) max(PT(1,:))];
%      minrow          maxrow      mincol       maxcol  

end

%%

% HOMOTRANS - homogeneous transformation of points
%
% Function to perform a transformation on homogeneous points/lines
% The resulting points are normalised to have a homogeneous scale of 1
%
% Usage:
%           t = homoTrans(P,v);
%
% Arguments:
%           P  - 3 x 3 or 4 x 4 transformation matrix
%           v  - 3 x n or 4 x n matrix of points/lines

%  Peter Kovesi
%  School of Computer Science & Software Engineering
%  The University of Western Australia
%  pk @ csse uwa edu au
%  http://www.csse.uwa.edu.au/~pk
%
%  April 2000
%  September 2007

function t = homoTrans(P,v)
    
    [dim,npts] = size(v);
    
    if ~all(size(P)==dim)
	error('Transformation matrix and point dimensions do not match');
    end

    t = P*v;  % Transform

    for r = 1:dim-1     %  Now normalise    
	t(r,:) = t(r,:)./t(end,:);
    end
    
    t(end,:) = ones(1,npts);
    
    
end

%%

function [ BB ] = calcbb( H, im )
%CALCBB Compute bounding box

[imwidth imheight tempx] = size(im);

% Bounding box calculation

xmax=0;
ymax=0;
xmin=0;
ymin=0;

% Transform Point(1,1)
	oldPoint = [1; 1; 1];
	newPoint= H* oldPoint;
	newPoint= newPoint/newPoint(3);

if (newPoint(1) < xmin) 
        xmin = newPoint(1);
end
if (newPoint(2) < ymin) 
        ymin = newPoint(2);
end
if (newPoint(1) > xmax) 
        xmax = newPoint(1);
end
if (newPoint(2) > ymax) 
        ymax = newPoint(2);
end
    
% Transform Point(1,imheight)
	oldPoint = [1; imheight; 1];
	newPoint= H* oldPoint;
	newPoint= newPoint/newPoint(3);

if (newPoint(1) < xmin) 
        xmin = newPoint(1);
end
if (newPoint(2) < ymin) 
        ymin = newPoint(2);
end
if (newPoint(1) > xmax) 
        xmax = newPoint(1);
end
if (newPoint(2) > ymax) 
        ymax = newPoint(2);
end

% Transform Point(imwidth,1)
	oldPoint = [imwidth; 1; 1];
	newPoint= H* oldPoint;
	newPoint= newPoint/newPoint(3);

if (newPoint(1) < xmin) 
        xmin = newPoint(1);
end
if (newPoint(2) < ymin) 
        ymin = newPoint(2);
end
if (newPoint(1) > xmax) 
        xmax = newPoint(1);
end
if (newPoint(2) > ymax) 
        ymax = newPoint(2);
end

% Transform Point(imwidth,imheight)
	oldPoint = [imwidth; imheight; 1];
	newPoint= H* oldPoint;
	newPoint= newPoint/newPoint(3);

if (newPoint(1) < xmin) 
        xmin = newPoint(1);
end
if (newPoint(2) < ymin) 
        ymin = newPoint(2);
end
if (newPoint(1) > xmax) 
        xmax = newPoint(1);
end
if (newPoint(2) > ymax) 
        ymax = newPoint(2);
end

BB = [xmin, xmax; ymin, ymax];

% End of BB calculation


end


