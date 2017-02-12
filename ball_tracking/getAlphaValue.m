function [value] = getAlphaValue(ellipse_t , x0 , y0 , img_grad_angle)
    a = ellipse_t.a;
    b = ellipse_t.b;
    phi = ellipse_t.phi;
    h = ellipse_t.X0_in; % changed the x,y properties, because that's how it appears in the row/col of image.
    k = ellipse_t.Y0_in;
    
    alpha = b ^ 2 * cos(phi) ^ 2 + a ^2 * sin(phi) ^2;
    beta = b ^ 2 * sin(phi) ^ 2 + a ^2 * cos(phi) ^2;
    
    gamma = sin(phi) * cos(phi) * (b ^2 - a^2);
    
    A = alpha;
    B = beta;
    C = -2*h*alpha -2*k*gamma;
    D = -2*k*beta - 2*h*gamma;
    E = 2*gamma;
    F = 2*h*k*gamma + alpha*h^2 + beta*k^2 - a^2 * b^2;
    
    % get the point on the ellipse from x0, and see which point is close to
    % y0.
    a_dash = B;
    b_dash = D + E*x0;
    c_dash = F + A*x0^2 + C*x0;
    
    discriminant = b_dash^2 - 4*a_dash*c_dash;
    if discriminant >=0
        y_one = (-b_dash + sqrt(discriminant))/2*a_dash;
        y_two = (-b_dash - sqrt(discriminant))/2*a_dash;
        if abs(y0 - y_one) > abs(y0 - y_two)
            y_point = y_two;
        else
            y_point = y_one;
        end
    else
        value = -1;
        return;
    end
    
    
    m_tangent = -1 * (2*A*x0 + C + E*y_point) / (2*B*y_point + D + E*x0);
    
    m_normal = -1 / m_tangent;
    
    normal_angle = atan(m_normal);
    
    value = abs(normal_angle - img_grad_angle);
    if value > pi
        value = value - pi;
    end
end