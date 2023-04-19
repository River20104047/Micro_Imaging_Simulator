function [centers,f1] = Fun01_PrintParticle_v01(Count, W, H, da, rgb)
%%
centers = zeros(Count, 2);

% generate the centers of the circles randomly inside the square
while true
    % generate random centers
    centers = rand(Count, 2) .* [W - da, H - da] + da/2;
    
    % check if all circles are inside the square
    if all(min(centers - da/2) > 0) && all(max(centers + da/2) < [W, H])
        break;
    end
end

% check if any of the two circles are overlapping with each other
overlapping = true;
while overlapping
    overlapping = false;
    for i = 1:Count
        for j = i+1:Count
            % calculate the distance between the centers of the circles
            d = norm(centers(i,:) - centers(j,:));
            % check if the circles are overlapping
            if d < da(i)/2 + da(j)/2
                % if overlapping, randomly reposition one of the circles
                centers(j,:) = rand(1,2) .* [W - da(j), H - da(j)] + da(j)/2;
                overlapping = true;
            end
        end
    end
end

% plot the circles

f1 = figure;

% % black background
% rectangle('Position', [0 - 0.1*W, 0 - 0.1*H, W + 0.2*W,H + 0.2*H], ...
%           'FaceColor', 'none','EdgeColor','black');
rectangle('Position', [0, 0, W,H], ...
          'FaceColor', 'none','EdgeColor',[0.7 0.7 0.7]);
hold on;
for i = 1:Count
    rectangle('Position', [centers(i,1)-da(i)/2, centers(i,2)-da(i)/2, da(i), da(i)], ...
        'Curvature', [1 1], 'FaceColor', rgb(i,:)/255,'EdgeColor','none');
end

axis equal;
xlim([0 - 0.1*W, W + 0.1*W]);
ylim([0 - 0.1*H, W + 0.1*H]);
xlabel('X [μm]')
ylabel('Y [μm]')
title('Distribution of particles')

end