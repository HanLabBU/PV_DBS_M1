% Store the speed
temp_speed = cur_speed;

% Find the artefact points
art_idx = find((temp_speed - mean(temp_speed))./std(temp_speed) > 6);

% Find the flanking points
x = [art_idx - 1; art_idx + 1];
x = x(:)';
y = temp_speed(x);

% Calculate interpolation points
qx = (x(1:2:end - 1) + x(2:2:end) )./ 2;

% Check that the values are the same
if art_idx == qx
    disp('Artefact and question interpolation points are not the same!!');
    temp_speed = [];
    return;
end

% calculate the individual points
vq = interp1(x, y, art_idx);

% Replace each point with the new points
temp_speed(art_idx) = vq;
