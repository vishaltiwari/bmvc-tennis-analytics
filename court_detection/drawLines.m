function [court_image] = drawLines(lines,width,height,frame_image)
%court_image = zeros(height,width);

%Try using the court_images along with the lines.
court_image = frame_image;

for k = 1:length(lines)
    xy = [lines(k).point1; lines(k).point2];
    x1 = xy(1,1);
    x2 = xy(2,1);
    
    y1 = xy(1,2);
    y2 = xy(2,2);

    m = Inf;
    c = Inf;

    if x1~=x2
        m = (y2-y1)/(x2-x1);
        c = y1 - (((y2-y1)*x1) / (x2-x1));
    end
    
    final_p = [];
    if m~=Inf
        %intersection with x=1
        ya = m * 1 + c; % xa = 1
        yb = m * width + c; % xb = width
        
        xc = (1 - c)/m; % yc = 1
        xd = (height - c)/m; %yd = height;
        
        pos = [1 ya ; width yb ; xc 1 ; xd height];
        
        for j=1:size(pos,1)
            p = pos(j,:);
            if p(1) >=1 && p(1) <=width && p(2) >=1 && p(2) <=height
                final_p = [final_p ; p];
            end
        end
    else
        %x = x0;
        final_p = [x1 1 ; x1 height];
    end
    
    %Draw the lines.
    court_image = insertShape(court_image , 'Line',[final_p(1,1) final_p(1,2) final_p(2,1) final_p(2,2)],'Color','white','LineWidth',3);
    %imshow(court_image);
end

end