function [] = read_and_filter_image()
I = imread("Apples.bmp");
G = mat2gray(I, [0 255]);
%21, 31, 51
plot_M_point_average_frequency_response_magnitude(21);
figure
imshow(blur_horizontal(G, 21));
plot_M_point_average_frequency_response_magnitude(31);
figure
imshow(blur_horizontal(G, 31));
plot_M_point_average_frequency_response_magnitude(51);
figure
imshow(blur_horizontal(G, 51));
plot_first_difference_frequency_response_magnitude();
figure;
imshow(highlight_vertical_lines(G));
end
function [] = plot_M_point_average_frequency_response_magnitude(L)
range = -pi:(pi/1000):pi;
bb = ones(1, L)/L;
HH = freqz(bb, 1, range);
figure
plot(range, abs(HH));
xlabel("w");
ylabel("|H(e^(jw))|");
title("Magnitude of frequency response for the " + L + " point average function");
end
%A is a 2d matrix
function filtered = blur_horizontal(A, m)
filtered = movmean(A, m, 2);
end
%A is a 2d matrix
function filtered = blur_vertical(A, m)
filtered = movmean(A, m, 1);
end
function [] = plot_first_difference_frequency_response_magnitude()
start_end = pi;
range = -start_end:(start_end/2000):start_end;
HH = 2 * abs(sin(range/2));
figure
plot(range, abs(HH));
xlabel("w");
ylabel("|H(e^(jw))|");
title("Magnitude of frequency response for the first difference function");
end
function filtered = highlight_vertical_lines(A)
filtered = filter([1, -1], 1, A, [], 2);%, dim=2);
end
function filtered = highlight_horizontal_lines(A)
filt
e
r
e
d
=
filt
e
r
(
[
1, -
1
], 1, A, [
], 1
);%, dim
=
2
); en