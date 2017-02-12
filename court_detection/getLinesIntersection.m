function [cols,rows] = getLinesIntersection(lines,width, height)
syms x y
eq_lst = [];
slopes = [];
cs = [];
for k = 1:length(lines)
    xy = [lines(k).point1; lines(k).point2];
    x1 = xy(1,1);
    x2 = xy(2,1);
    
    y1 = xy(1,2);
    y2 = xy(2,2);

    m = Inf;
    c = Inf;
    if x1 == x2
        eq = y == x1;
    end
    if x1~=x2
        eq = y - ((y2 - y1)/(x2-x1))*x == y1 - (((y2-y1)*x1) / (x2-x1));
        m = (y2-y1)/(x2-x1);
        c = y1 - (((y2-y1)*x1) / (x2-x1));
    end
    slopes = [slopes m];
    cs = [cs c];
    eq_lst = [eq_lst eq];
    
end

% Solve each pair
cols = [];
rows = [];
%hold on;
for i=1:size(eq_lst,2)
    eq1 = eq_lst(i);
%     xy = [lines(i).point1; lines(i).point2];
%     plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');
    for j=i+1:size(eq_lst,2)
        if slopes(i) == slopes(j)
            continue;
        end
                
%         xy2 = [lines(j).point1; lines(j).point2];
%         plot(xy2(:,1),xy2(:,2),'LineWidth',2,'Color','green');
        
        eq2 = eq_lst(j);
        soln = solve([eq1 , eq2] , [x,y]);
        col = soln.x;
        row = soln.y;
        if isempty(col) || isempty(row)
            continue;
        end
        if col>=1 && col <=width && row>=1 && row<=height
            cols = [cols ; double(col)];
            rows = [rows ; double(row)];
        end
    end
end

end