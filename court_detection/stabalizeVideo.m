function [] = stabalizeVideo(video_in , video_out)
    % Reset the video source to the beginning of the file.
    
    %hVideoSrc = VideoReader(video_in);
    %reset(hVideoSrc);

    hVideoSrc = vision.VideoFileReader(video_in);
    % Process all frames in the video
    movMean = step(hVideoSrc);
    imgB = movMean;
    imgBp = imgB;
    correctedMean = imgBp;
    ii = 2;
    Hcumulative = eye(3);
    while ~isDone(hVideoSrc) && ii < 10
        % Read in new frame
        imgA = imgB; % z^-1
        imgAp = imgBp; % z^-1
        imgB = step(hVideoSrc);
        movMean = movMean + imgB;

        % Estimate transform from frame A to frame B, and fit as an s-R-t
        H = cvexEstStabilizationTform(imgA,imgB);
        HsRt = cvexTformToSRT(H);
        Hcumulative = HsRt * Hcumulative;
        imgBp = imwarp(imgB,affine2d(Hcumulative),'OutputView',imref2d(size(imgB)));

        % Display as color composite with last corrected frame
        step(hVPlayer, imfuse(imgAp,imgBp,'ColorChannels','red-cyan'));
        correctedMean = correctedMean + imgBp;

        ii = ii+1;
    end
    correctedMean = correctedMean/(ii-2);
    movMean = movMean/(ii-2);

    % Here you call the release method on the objects to close any open files
    % and release memory.
    release(hVideoSrc);
    release(hVPlayer);
end