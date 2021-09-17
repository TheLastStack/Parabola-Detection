clear all;
img = imread('parabola.bmp');
% The read image interprets white pixels as 255 and black pixels as 0
[yl, xl] = size(img);
% Getting dimensions of the image
xlength = 301;
ylength = 301;
% Dimensions of vertex search space
xcenter = 4;
ycenter = 5;
% Center point of vertex search space
xvx = 1:xlength;
% Defining first dimension of accumulator cells: vertex x-position
yvy = 1:ylength;
% Defining second dimension of accumulator cells: vertex y-position
p_acc = zeros(ylength, xlength, floor(sqrt(xlength^2 + ylength^2)) - 1);
% Defining third dimension of accumulator cells: Distance between focus and
% vertex. The maximum detectable value this dimension can take is sqrt(xlength**2 +
% ylength**2) - 1.

for i = 1:xl
    for j = 1:yl
        %Looping over all points in image
        if img(j, i) ~= 0
            %For thos points which are white, we graph it in accumulator
            %cells
            for iv = xvx
                for jv = yvy
                    %Graphing over first and second dimensions. Together
                    %they represent the location of vertex
                    if jv ~= j
                        opx = iv - (floor(xlength / 2) + 1 - xcenter);
                        opy = jv - (floor(ylength / 2) + 1 - ycenter);
                        % Potential vector coordinates in image coordinate
                        % system
                        opz = round(((i - opx) ^ 2)/(4 * (j - opy)));
                        if opz > 0 && opz <= floor(sqrt(xlength^2 + ylength^2)) - 1
                            p_acc(jv, iv, opz) = p_acc(jv, iv, opz) + 1;
                        end
                    end
                end
            end
        end
    end
end
[max_value, max_index] = max(p_acc(:));
all_max_indices = find(p_acc == max_value);
[vertex_y, vertex_x, p] = ind2sub(size(p_acc), round(mean(all_max_indices)));
% In the above steps, we assume there is only one parabola in the image.
% We search for the best vertex position by averaging all of the best
% possible vertex positions
vertex_y = vertex_y - (floor(ylength / 2) + 1 - ycenter);
vertex_x = vertex_x - (floor(xlength / 2) + 1 - xcenter);
% This gives vertex position in the coordinate system of the image 
imgcopy = img;
% Copying the image
display_xlength = 601;
display_ylength = 601;
display_xcenter = 0;
display_ycenter = 0;
display_vertex_y = vertex_y - display_ycenter;
display_vertex_x = vertex_x - display_xcenter;
plot_points_x = -floor(display_xlength/2) - display_xcenter:floor(display_xlength/2) - display_xcenter;
plot_points_y = -floor(display_ylength/2) - display_ycenter:floor(display_ylength/2) - display_ycenter;
plot_parabola = zeros(1, length(plot_points_x));
for i = 1:display_xlength
    j = ((plot_points_x(i) - display_vertex_x)^2/(4 * p)) + display_vertex_y;
    if j < max(plot_points_y) && j > min(plot_points_y)
        plot_parabola(i) = j;
    end
end
figure;
plot(plot_points_x, plot_parabola);
hold on;
image_parabola = zeros(1, length(plot_points_x));
for i = 1:display_xlength
    if plot_points_x(i) >= 1 && plot_points_x(i) <= xl
        [val,j] = max(img(:, plot_points_x(i)) > 0);
        if val > 0 && j > min(plot_points_y) && j < max(plot_points_y)
            image_parabola(i) = j;
        end
    end
end
plot(plot_points_x, image_parabola);
% Overlapping parabola in image with parabola obtained by plotting with
% obtained vertex co-ordinates and distance between focus and vertex (p)
% Where values are not available/ values go outside the display window, the
% plot takes on 0 value.
% Please note that in the image co-ordinate system, y increases downwards,
% so the parabola opens in the direction of increasing y. So, when plotted
% with increasing y pointing upwards, the shape of the parabola is
% preserved, but the image 'flips'.
meter_to_pixel_ratio = 5;
% Scaling of image
% Assuming it to be 5m for each pixel here.
%xo = <xval>;
%yo = <yval>;
% Give starting point of trajectory.
xo = 1;
[~, yo] = max(img(:,1)>0);
yo = -yo + yl; %(Inverting and shifting y axis)
calc_x = vertex_x - 1;
xo = xo - 1;
% Shifting x-axis to x=1
% Here, we assume the projectile starts at the left most point in the image
% This is a default value and can be replaced
% Please note that this is in pixel coordinates
theta = atan(calc_x / (2 * p) + 10 * xo / (sqrt(20 * p)));
initial_velocity = sqrt( 20 * p * meter_to_pixel_ratio) / cos(theta);
% Returns velocity and angle
theta * 180 / pi



        




            
        
