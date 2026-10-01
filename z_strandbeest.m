clear()
set(groot, 'defaultTextInterpreter', 'latex')
set(groot, 'defaultAxesTickLabelInterpreter', 'latex')
set(groot, 'defaultLegendInterpreter', 'latex')

leg_params = define_leg_parameters();
solver_params = struct();


%%% run animation

theta = 0;
vertex_coords_guess = [0,50,-50,0,-50,50,-100,0,-100,-50,-50,-50,-50,-100];
leg_params.leg_pos = compute_coords(vertex_coords_guess, leg_params,theta);
leg_params.leg_velocity = compute_velocities(vertex_coords_guess,leg_params,theta);
leg_drawing = initialize_leg_drawing(leg_params);

complete_vertex_coords = zeros(14,2500);
complete_leg_velocity = zeros(14,2);
complete_leg_velocity_jacobian = zeros(14,2);
theta = linspace(0,10*pi,length(complete_vertex_coords));
for i = 1:length(complete_vertex_coords)
    complete_vertex_coords(:,i)=compute_coords(vertex_coords_guess,leg_params,theta(i));
    [complete_leg_velocity(i,:),complete_leg_velocity_jacobian(i,:)]=compute_velocities(vertex_coords_guess,leg_params,theta(i));
end
update_leg_drawing(complete_vertex_coords,leg_drawing,leg_params,complete_leg_velocity);
plot_velocities(theta, complete_leg_velocity, complete_leg_velocity_jacobian)


%Defines a struct that contains the constants for Strandbeest leg
%OUTPUTS:
%leg_params: a struct with the following fields
%leg_params.num_vertices: number of vertices in linkage
%leg_params.num_linkages: number of links in linkage
%leg_params.link_lengths: list of lengths for each link in the leg mechanism
%leg_params.crank_length: length of crank shaft
%leg_params.vertex_pos0: fixed position coords of vertex 0
%leg_params.vertex_pos2: fixed position coords of vertex 2
function leg_params = define_leg_parameters()
%initialize leg_params structure
leg_params = struct();
%number of vertices in linkage
leg_params.num_vertices = 7;
%number of links in linkage
leg_params.num_linkages = 10;
%matrix relating links to vertices
leg_params.link_to_vertex_list = ...
[ 1, 3;... %link 1 adjacency
3, 4;... %link 2 adjacency
2, 3;... %link 3 adjacency
2, 4;... %link 4 adjacency
4, 5;... %link 5 adjacency
2, 6;... %link 6 adjacency
1, 6;... %link 7 adjacency
5, 6;... %link 8 adjacency
5, 7;... %link 9 adjacency
6, 7 ... %link 10 adjacency
];
%list of lengths for each link
%in the leg mechanism
leg_params.link_lengths = ...
[ 50.0,... %link 1 length
55.8,... %link 2 length
41.5,... %link 3 length
40.1,... %link 4 length
39.4,... %link 5 length
39.3,... %link 6 length
61.9,... %link 7 length
36.7,... %link 8 length
65.7,... %link 9 length
49.0 ... %link 10 length
];
%length of crank shaft
leg_params.crank_length = 15.0;
%fixed position coords of vertex 0
leg_params.vertex_pos0 = [0;0];
%fixed position coords of vertex 2
leg_params.vertex_pos2 = [-38.0;-7.8];
end

%% Computes stuff
%Computes the vertex coordinates that describe a legal linkage configuration
%INPUTS:
%vertex_coords_guess: a column vector containing the (x,y) coordinates of every vertex
% these coords are just a GUESS! It's used to seed Newton's method
%leg_params: a struct containing the parameters that describe the linkage
%theta: the desired angle of the crank
%OUTPUTS:
%vertex_coords_root: a column vector containing the (x,y) coordinates of every vertex
% these coords satisfy all the kinematic constraints!
function vertex_coords_root = compute_coords(vertex_coords_guess, leg_params,theta)
%your code here
%you will likely need to make a wrapper function of linkage_error_func
%so that it is only a function of vertex_coords
%(and not leg_params or theta, which should be set beforehand)
%you can then pass this wrapper function to your multidimensional Newton
%solver, along with vertex_coords_guess to find the vertex coordinates
%corresponding to the legal configuration of the linkage,
%given the values set for leg_params and theta
solver_params = struct();
vertex_coords_root = multi_newton_solver(@(vertex_coords) linkage_error_func(vertex_coords,leg_params,theta),vertex_coords_guess, solver_params);
end

%Computes the theta derivatives of each vertex coordinate for the Jansen linkage
%INPUTS:
%vertex_coords: a column vector containing the (x,y) coordinates of every vertex
% these are assumed to be legal values that are roots of the error funcs!
%leg_params: a struct containing the parameters that describe the linkage
%theta: the current angle of the crank
%OUTPUTS:
%dVdtheta: a column vector containing the theta derivates of each vertex coord
function [dVdtheta_explicit,dVdtheta_jacobian] = compute_velocities(vertex_coords, leg_params, theta)
    cvc = @(theta) compute_coords(vertex_coords, leg_params, theta);
    dtheta = 1e-7;
    ddtheta_positive = column_to_matrix(cvc(theta+dtheta));
    ddtheta_negative = column_to_matrix(cvc(theta-dtheta));
    dxdtheta_explicit = (ddtheta_positive(7,1)-ddtheta_negative(7,1))/(2*dtheta);
    dydtheta_explicit = (ddtheta_positive(7,2)-ddtheta_negative(7,2))/(2*dtheta);
    dVdtheta_explicit = [dxdtheta_explicit,dydtheta_explicit];
    
    vertex_coords_true = cvc(theta);
    J = approximate_jacobian(@(vertex_coords) link_length_error_func(vertex_coords,leg_params),vertex_coords_true);
    M = [eye(4),zeros(4,10);J];
    B = [-leg_params.crank_length*sin(theta);leg_params.crank_length*cos(theta);zeros(12,1)];
    
    dvdtheta = column_to_matrix(M\B);
    dVdtheta_jacobian = [dvdtheta(7,1),dvdtheta(7,2)];

end


%% Error functions
%Error function that encodes all necessary linkage constraints
%INPUTS:
%vertex_coords: a column vector containing the (x,y) coordinates of every vertex
%leg_params: a struct containing the parameters that describe the linkage
%theta: the current angle of the crank
%OUTPUTS:
%error_vec: a vector describing each constraint on the linkage
% when error_vec is all zeros, the constraints are satisfied
function error_vec = linkage_error_func(vertex_coords, leg_params, theta)
    distance_errors = link_length_error_func(vertex_coords, leg_params);
    coord_errors = fixed_coord_error_func(vertex_coords, leg_params, theta);
    error_vec = [distance_errors;coord_errors];
end


%Error function that encodes the link length constraints
%INPUTS:
%vertex_coords: a column vector containing the (x,y) coordinates of every vertex
% in the linkage. There are two ways that I would recommend stacking
% the coordinates. You could alternate between x and y coordinates:
% i.e. vertex_coords = [x1;y1;x2;y2;...;xn;y_n], or alternatively
% you could do all the x's first followed by all of the y's
% i.e. vertex_coords = [x1;x2;...xn;y1;y2;...;yn]. You could also do
% something else entirely, the choice is up to you.
%leg_params: a struct containing the parameters that describe the linkage
% importantly, leg_params.link_lengths is a list of linakge lengths
% and leg_params.link_to_vertex_list is a two column matrix where
% leg_params.link_to_vertex_list(i,1) and
% leg_params.link_to_vertex_list(i,2) are the pair of vertices connected
% by the ith link in the mechanism
%OUTPUTS:
%length_errors: a column vector describing the current distance error of the ith
% link specifically, length_errors(i) = (xb-xa)^2 + (yb-ya)^2 - d_i^2
% where (xa,ya) and (xb,yb) are the coordinates of the vertices that
% are connected by the ith link, and d_i is the length of the ith link
function length_errors = link_length_error_func(vertex_coords, leg_params)
    length_errors = zeros(leg_params.num_linkages,1);
    % Debugging vertex_coords format
    % vertex_coords
    vertex_coords_matrix = column_to_matrix(vertex_coords);
    for i=1:leg_params.num_linkages
        current_vertex = leg_params.link_to_vertex_list(i,:);
        xb = vertex_coords_matrix(current_vertex(2),1); % gets the x coordinate of the second vertex
        xa = vertex_coords_matrix(current_vertex(1),1); % gets the x coordintae of the first vertex
        yb = vertex_coords_matrix(current_vertex(2),2); % gets the y coordinate of the second vertex
        ya = vertex_coords_matrix(current_vertex(1),2); % gets the y coordinate of the first vertex
        length_errors(i) = sqrt((xb-xa)^2+(yb-ya)^2)-leg_params.link_lengths(i);
    end
end


%Error function that encodes the fixed vertex constraints
%INPUTS:
%vertex_coords: a column vector containing the (x,y) coordinates of every vertex
% same input as link_length_error_func
%leg_params: a struct containing the parameters that describe the linkage
% importantly, leg_params.crank_length is the length of the crank
% and leg_params.vertex_pos0 and leg_params.vertex_pos2 are the
% fixed positions of the crank rotation center and vertex 2.
%theta: the current angle of the crank
%OUTPUTS:
%coord_errors: a column vector of height four corresponding to the differences
% between the current values of (x1,y1),(x2,y2) and
% the fixed values that they should be
function coord_errors = fixed_coord_error_func(vertex_coords, leg_params, theta)
    x1 = vertex_coords(1); y1 = vertex_coords(2); 
    x2 = vertex_coords(3); y2 = vertex_coords(4);

    bar_x1 = leg_params.vertex_pos0(1) + leg_params.crank_length*cos(theta);
    bar_y1 = leg_params.vertex_pos0(2) + leg_params.crank_length*sin(theta);
    bar_x2 = leg_params.vertex_pos2(1);
    bar_y2 = leg_params.vertex_pos2(2);
    coord_errors = [x1-bar_x1,...
                    y1-bar_y1,...
                    x2-bar_x2,...
                    y2-bar_y2].';
end

%% Matrix conversion

%Converts from the column vector form of the coordinates to a
%friendlier matrix form
%INPUTS:
%coords_in = [x1;y1;x2;y2;...;xn;yn] (2n x 1 column vector)
%OUTPUTS:
%coords_out = [x1,y1;x2,y2;...;xn,yn] (n x 2 matrix)
function coords_out = column_to_matrix(coords_in)
num_coords = length(coords_in);
coords_out = [coords_in(1:2:(num_coords-1)),coords_in(2:2:num_coords)];
end
%Converts from the matrix form of the coordinates back to the
%original column vector form
%INPUTS:
%coords_in = [x1,y1;x2,y2;...;xn,yn] (n x 2 matrix)
%OUTPUTS:
%coords_out = [x1;y1;x2;y2;...;xn;yn] (2n x 1 column vector)
function coords_out = matrix_to_column(coords_in)
num_coords = 2*size(coords_in,1);
coords_out = zeros(num_coords,1);
coords_out(1:2:(num_coords-1)) = coords_in(:,1);
coords_out(2:2:num_coords) = coords_in(:,2);
end



%% Drawing
%Creates a set of plotting objects to keep track of each link drawing
%each vertex drawing, and the crank drawing
%INPUTS:
%leg_params: a struct containing the parameters that describe the linkage
%OUTPUTS:
%leg_drawing: a struct containing all the plotting objects for the linkage
% leg_drawing.linkages is a cell array, where each element corresponds
% to a plot of a single link (excluding the crank)
% leg_drawing.crank is a plot of the crank link
% leg_drawing.vertices is a cell array, where each element corresponds
% to a plot of one of the vertices in the linkage
function leg_drawing = initialize_leg_drawing(leg_params)
leg_drawing = struct();
leg_drawing.fig1 = figure('units','pixels','position',[0 0 1440 1056]);hold on;
leg_drawing.linkages = cell(leg_params.num_linkages,1);
leg_pos=column_to_matrix(leg_params.leg_pos);
for linkage_index = 1:leg_params.num_linkages
    leg_index = leg_params.link_to_vertex_list(linkage_index,:);
    leg_drawing.linkages{linkage_index} = line([leg_pos(leg_index,1)],[leg_pos(leg_index,2)],'color','k','linewidth',2,'HandleVisibility','off');
end
leg_drawing.crank = line([0,0],[0,0],'color','magenta','linewidth',1.5,'HandleVisibility','off');
leg_drawing.vertices = cell(leg_params.num_vertices,1);
for vertex_index = 1:leg_params.num_vertices
    leg_drawing.vertices{vertex_index} = line(leg_pos(vertex_index,1),leg_pos(vertex_index,2),'marker',...
'o','markerfacecolor','r','markeredgecolor','r','markersize',8,'HandleVisibility','off',LineStyle='none');
end
leg_drawing.foot_trace = line([0,0],[0,0],'linewidth',1.5,'HandleVisibility','off',Color="blue",LineStyle=":");
[x_coords, y_coords, xhead, yhead] = directionLine(leg_pos(7,1),leg_pos(7,2),leg_params.leg_velocity(1),leg_params.leg_velocity(2),.8);
leg_drawing.leg_velocity = line([x_coords,xhead],[y_coords,yhead],'linewidth',1.5,'HandleVisibility','off',Color="red");
end




function update_leg_drawing(complete_vertex_coords, leg_drawing, leg_params, complete_leg_velocity, record)
arguments
    complete_vertex_coords
    leg_drawing
    leg_params
    complete_leg_velocity = 0
    record = false
end
if record
    mypath = 'C:\Documents\meche_math\assignment_2\';
    fname = 'strandbeest_animation.avi';
    input_fname = [mypath,fname];

    %create a videowriter, which will write frames to the animation file
    writerObj = VideoWriter(input_fname);
    open(writerObj); %must call open before writing any frames
end
leg_drawing.fig1
title('Strandbeest Animation');
ylabel('$y (-)$');
xlabel('$x (-)$');
fontsize(scale=2);
axis([-120,20,-100,40]); 
line([0,complete_vertex_coords(3,1)],[0,complete_vertex_coords(4,1)],'marker',...
'o','markerfacecolor','g','markeredgecolor','g','markersize',8,'HandleVisibility','off',LineStyle='none');
%generate legend
line(0,0,'color','k','linewidth',2);
line(0,0,'marker','o','markerfacecolor','r','markeredgecolor','r','markersize',8,LineStyle='none')
line(0,0,'color','magenta','linewidth',1.5);
line(0,0,'linewidth',1.5,Color="blue",LineStyle=":")
line(0,0,'marker','o','markerfacecolor','g','markeredgecolor','g','markersize',8,LineStyle='none')
line(0,0,Color="red");
legend('Linkages','Vertices','Crank','Foot Path','Fixed Points','Leg Velocity Overlay','Location','northeastoutside');



%complete_vertex_coords is a [14 x t] matrix
for i = 1:length(complete_vertex_coords)
    vertex_coord = complete_vertex_coords(:,i);
    vertex_coord = column_to_matrix(vertex_coord);

    %iterate through each link, and update corresponding link plot
    for linkage_index = 1:leg_params.num_linkages
    %linkage_index is the label of the current link
    %your code here
    %line_x and line_y should both be two element arrays containing
    %the x and y coordinates of the line segment describing the current link
        leg_index = leg_params.link_to_vertex_list(linkage_index,:);
        line_x = [vertex_coord(leg_index,1)];
        line_y = [vertex_coord(leg_index,2)];
        set(leg_drawing.linkages{linkage_index},'xdata',line_x,'ydata',line_y);
    end
    %iterate through each vertex, and update corresponding vertex plot


    for vertex_index = 1:leg_params.num_vertices
        line_x = vertex_coord(:,1);
        line_y = vertex_coord(:,2);
        set(leg_drawing.vertices{vertex_index},'xdata',line_x,'ydata',line_y);
        
    end
    %your code here
    %crank_x and crank_y should both be two element arrays
    %containing the x and y coordinates of the line segment describing the crank
    crank_x = [0,complete_vertex_coords(1,i)];
    crank_y = [0,complete_vertex_coords(2,i)];
    set(leg_drawing.crank,'xdata',crank_x,'ydata',crank_y);
    set(leg_drawing.foot_trace,'xdata',complete_vertex_coords(end-1,1:i),'ydata',complete_vertex_coords(end,1:i));
    if complete_leg_velocity
        
        [x_coords, y_coords, xhead, yhead] = directionLine(vertex_coord(7,1),vertex_coord(7,2),complete_leg_velocity(i,1),complete_leg_velocity(i,2),.8);
        set(leg_drawing.leg_velocity,'xdata',[x_coords,xhead],'ydata',[y_coords,yhead]);

        % set(leg_drawing.leg_velocity,'xdata',[vertex_coord(7,1),vertex_coord(7,1)+.5*complete_leg_velocity(i,1)],'ydata',[vertex_coord(7,2),vertex_coord(7,2)+.5*complete_leg_velocity(i,2)]);
        
    end
    
    drawnow;
    if record
        %capture a frame (what is currently plotted)
        current_frame = getframe(leg_drawing.fig1);
        %write the frame to the video
        writeVideo(writerObj,current_frame);
    end
end
if record
    current_frame = getframe(leg_drawing.fig1);
    %write the frame to the video
    writeVideo(writerObj,current_frame);
    close(writerObj);
end
end

function [x_coords, y_coords, x_head, y_head] = directionLine(x0, y0, x_direction, y_direction, scale_factor, head_length, head_angle_deg)
% DIRECTIONLINE Generate line + arrowhead coordinates from a point along a direction
%
% The endpoint is computed as (x0,y0) + scale_factor * (x_direction, y_direction),
% so scale_factor is a raw multiplier on the direction vector (no normalization).
%
% Inputs:
%   x0, y0           - Starting point coordinates (scalars)
%   x_direction      - X component of direction vector
%   y_direction      - Y component of direction vector
%   scale_factor     - Multiplier applied to the direction vector
%   head_length      - (optional) Length of arrowhead barbs. Default: 10% of shaft length
%   head_angle_deg   - (optional) Half-angle of arrowhead in degrees. Default: 20
%
% Outputs:
%   x_coords, y_coords - Shaft coordinates for line(x_coords, y_coords)
%   x_head, y_head     - Arrowhead coordinates (3 points per barb, NaN-separated)
%                        ready for line(x_head, y_head)
%
% Example:
%   [xc, yc, xh, yh] = directionLine(0, 0, 1, 1, 5);
%   figure; axis equal; hold on;
%   line(xc, yc, 'Color', 'r', 'LineWidth', 2);
%   line(xh, yh, 'Color', 'r', 'LineWidth', 2);

    % --- Validate direction ---
    dir_mag = hypot(x_direction, y_direction);
    if dir_mag == 0
        error('Direction vector cannot be zero.');
    end

    % Unit direction (only used for orienting the arrowhead barbs,
    % NOT for scaling the shaft)
    ux = x_direction / dir_mag;
    uy = y_direction / dir_mag;

    % Shaft endpoints: raw multiplier on the direction vector
    x_start = x0;
    y_start = y0;
    x_end   = x0 + scale_factor * x_direction;
    y_end   = y0 + scale_factor * y_direction;

    x_coords = [x_start, x_end];
    y_coords = [y_start, y_end];

    % Actual shaft length (for default head sizing)
    shaft_length = hypot(x_end - x_start, y_end - y_start);

    % --- Defaults ---
    if nargin < 6 || isempty(head_length)
        head_length = 0.1 * shaft_length;
    end
    if nargin < 7 || isempty(head_angle_deg)
        head_angle_deg = 20;
    end

    % --- Arrowhead barbs ---
    theta = deg2rad(head_angle_deg);

    % Barb 1: rotate unit direction by +theta
    bx1 = x_end - head_length * (cos(theta)*ux - sin(theta)*uy);
    by1 = y_end - head_length * (sin(theta)*ux + cos(theta)*uy);

    % Barb 2: rotate unit direction by -theta
    bx2 = x_end - head_length * (cos(theta)*ux + sin(theta)*uy);
    by2 = y_end - head_length * (-sin(theta)*ux + cos(theta)*uy);

    % NaN-separated segments so a single line() call draws both barbs
    x_head = [x_end, bx1, x_end, NaN, x_end, bx2, x_end];
    y_head = [y_end, by1, y_end, NaN, y_end, by2, y_end];
end

function plot_velocities(theta, complete_leg_velocity,complete_leg_velocity_jacobian)
    figure();
    plot(theta,complete_leg_velocity(:,1),'b'); hold on;
    plot(theta,complete_leg_velocity_jacobian(:,1),'rx');
    title('Leg Velocity ($x$) over $\theta$','FontSize',36);
    xlabel('$\theta$ (radians)','FontSize',24);
    ylabel('$dx_{tip}/d\theta$ (-)','FontSize',24)
    legend('Finite Differences','Jacobian','FontSize',24,'Location','southeast')
    ax = gca;
    ax.FontSize=20;

    
    figure();
    plot(theta,complete_leg_velocity(:,2),'b'); hold on;
    plot(theta,complete_leg_velocity_jacobian(:,2),'rx');
    title('Leg Velocity ($y$) over $\theta$','FontSize',36);
    xlabel('$\theta$ (radians)','FontSize',24);
    ylabel('$dy_{tip}/d\theta$ (-)','FontSize',24)
    legend('Finite Differences','Jacobian','FontSize',24,'Location','southeast')
    ax = gca;
    ax.FontSize=20;
end