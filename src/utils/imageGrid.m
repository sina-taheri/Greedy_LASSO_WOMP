function Igrid = imageGrid(I, varargin)
% imageGrid puts a grid on an image
% Example 1: 
% I = imread('rice.png');
% imageGrid(I);
% Example 2:
% I = imread('cameraman.tif');
% Igrid = imageGrid(I, 'vertLines', 10, 'horzLines', 10, 'lineWidth', 1, 'lineStyle', ':', 'method', 'burn');
% figure
% imshowpair(I, Igrid, 'montage')
% Setting 'method' to 'draw' (default) opens a new figure with the grid
% drawn onto it
% Setting 'method' to 'burn' outputs a new image with the grid burned onto
% it, which requires the computer vision toolbox
% Tested in MATLAB R2018b
% Requires computer vision toolbox
% (c) Clay Swackhamer
% swackhamerclay at gmail.com
% Version 1.2 2018-10-14
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Set defaults
defaultVertLines = 10;
defaultHorzLines = 10;
defaultLineWidth = 1;
defaultLineStyle = '-';
defaultMethod = 'draw';
p = inputParser;
p.CaseSensitive = false;
p.KeepUnmatched = true; %store extra inputs instead of throwing error
addRequired(p, 'I');
addOptional(p, 'vertLines', defaultVertLines);
addOptional(p, 'horzLines', defaultHorzLines);
addOptional(p, 'lineWidth', defaultLineWidth);
addOptional(p, 'lineStyle', defaultLineStyle);
addOptional(p, 'method', defaultMethod);
%Parse inputs
parse(p, I, varargin{:});
vertLines = p.Results.vertLines;
horzLines = p.Results.horzLines;
lineWidth = p.Results.lineWidth;
lineStyle = p.Results.lineStyle;
method = p.Results.method;
drawLogical = strcmpi(method, 'draw'); %True if 'technique' is 'draw'
%Find out where to draw lines
horzTicks = round(linspace(1, size(I,1), horzLines),0); %cols where we will draw vertical lines
vertTicks = round(linspace(1, size(I,2), vertLines),0); %rows where we will draw horizontal lines
maxY = size(I,1); %num pixels in Y
maxX = size(I,2); %num pixels in X
if drawLogical
    %Display the image and draw lines on it
    imshow(I)
    for i = 1:1:vertLines
        line([vertTicks(i), vertTicks(i)], [maxY, 0], 'LineWidth', lineWidth, 'Color', 'k', 'LineStyle', lineStyle); %draw vertical lines
    end
    for i = 1:1:horzLines
        line([0, maxX], [horzTicks(i), horzTicks(i)], 'LineWidth', lineWidth, 'Color', 'k', 'LineStyle', lineStyle); %draw horizontal lines
    end
    Igrid = 'Lines drawn on open figure';
else
    
    %Burn the lines into the image (changes pixels)
    %Put together line coordinates the way that ShapeInserter wants 
    % [x1_line1, y1_line1, x2_line1, y2_line1;...
    %  x1_line2, y1_line2, x2_line2, y2_line2]
    %Vertical lines
    Vlines = zeros(length(vertTicks),4);
    Vlines(:,1) = vertTicks;
    Vlines(:,2) = 0;
    Vlines(:,3) = vertTicks;
    Vlines(:,4) = maxY;
    %Horizontal lines
    Hlines = zeros(length(horzTicks), 4);
    Hlines(:,1) = 0;
    Hlines(:,2) = horzTicks;
    Hlines(:,3) = maxX;
    Hlines(:,4) = horzTicks;
    %Put all lines together
    lines = vertcat(Vlines, Hlines);
    lines = int32(lines);
    %Insert lines and output the modified image
    Igrid = insertShape(I, 'Line', lines, 'LineWidth', lineWidth, 'Color', 'black');
end
end