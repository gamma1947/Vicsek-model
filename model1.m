%% Building the vicsek model

%% Initial setup
% Creating a lattice of random unit vectors

N = 5;
L = 1;
x = L.*rand(1,N);
y = L.*rand(1,N);
eta = 0.1;
velocity = 0.03;

scatter(x, y, "black.")
hold on
%% 
% Adding vectors to them-assigning random directions

theta =  2*pi.*rand(1,N);
u = cos(theta);
b = sin(theta);
quiver(x, y, u, b, 0, 'r'); % Plot unit vectors in red
hold off; % Release the plot


%% Get the wheel of time turning
% A function to update direction based on average direction of neighbours within 
% a radius
%% trial runs

neighbours_of_i = find_neighbours(5, x, y, 10, N);
[theta_average, updated_x, updated_y] = diksuchi(5, x, y, u, v, theta, neighbours_of_i, eta, velocity);

%% Simulation over time 

T = 1; % number of time steps
r = 1;
epsilon = 0.05;

x_t_i = x;
y_t_i = y;
theta_t_i = theta;

figure;
xlim([0 L]);
ylim([0 L]);

video_filename = 'vicsek.avi';
writer = VideoWriter(video_filename); 
open(writer);
metrics = zeros(1,T);

for j = 1:T
    for i = 1:N
        ngh = find_neighbours(i, x_t_i, y_t_i, r, N);
        [theta_average, updated_x, updated_y] = diksuchi(i, x_t_i, y_t_i, u, b, theta, ngh, eta, velocity);
        if updated_x > L
            x_t_i(i) = updated_x - L;
        elseif updated_x < epsilon
            x_t_i(i) = updated_x + L;
        else 
            x_t_i(i) = updated_x;
        end
        if updated_y > L
            y_t_i(i) = updated_y - L;
        elseif updated_x < epsilon
            y_t_i(i) = updated_y + L;
        else 
            y_t_i(i) = updated_y;
        end
        theta_t_i(i) = theta_average;
    end
    u_t_i = cos(theta_t_i);
    v_t_i = sin(theta_t_i);

    metrics(j) = metric(N, velocity, u_t_i, v_t_i);

    cla;
    xlim([0 L]);
    ylim([0 L]);
    scatter(x_t_i, y_t_i, "black.");
    hold on;
    quiver(x_t_i, y_t_i, u_t_i, v_t_i, 0, 'r');
    title(['Time step:' num2str(j)]);
    xlabel('X');
    ylabel('Y');
    drawnow;

    frame = getframe(gcf);
    writeVideo(writer, frame); % Use the new variable name 'writer'
end

close(writer); 
hold off

%%
t = 1:T ;
scatter(t, metrics)
hold off

%% Funtions
function neighbours = find_neighbours(poi_n, x, y, r, N)
    %poi = particle of interest
    x_n = zeros(1, N);
    y_n = zeros(1, N);
    poi_x = x(poi_n);
    poi_y = y(poi_n);
    n_coords = zeros(1, N);
    for i = 1:N
        if (poi_y-r)< y(i) && y(i) < (poi_y+r)
            if (poi_x-r)< x(i) && x(i) < (poi_x+r)
                dist = sqrt((y(i)-poi_y)^2 + (x(i)-poi_x)^2);
                if dist < r
                    y_n(i) = y(i);
                    x_n(i) = x(i);
                    n_coords(i) = i;
                end
            end
        end
    end
    neighbours = n_coords(n_coords ~= 0);
end

function [theta_avg, x_t1, y_t1] = diksuchi(poi_n, x, y, u, v,theta, neighbours, eta, velocity)
    n = neighbours;
    no_n = numel(n);
    u_n = zeros(1,no_n);
    v_n = zeros(1,no_n);
    for i = 1:no_n
        u_n(i) = u(n(i));
        v_n(i) = v(n(i));
    end
    noise = -eta/2 + eta*rand(1);
    a = mean(v_n);
    b = mean(u_n);
    theta_avg = atan(a/b) + noise;
    x_t1 = x(poi_n) + velocity*cos(theta(poi_n))*1;
    y_t1 = y(poi_n) + velocity*sin(theta(poi_n))*1;
end

function va = metric(N, velocity, u, v)
    va = (1/(N*velocity))*(abs(sum(u)+ sum(v)));
end