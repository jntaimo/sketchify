clear all; close all;
image = "img/imroom.jpg";
RGB = imread(image); % import image
figure(1); imshow(RGB); title("Original Image");
sharp = imsharpen(RGB, "radius", 10); % can adjust sharpen
gray = rgb2gray(sharp); % puts the image in grayscale

% "Canny" works better for more realistic / detailed pictures
% use default "Sobel" otherwise
edges = edge(gray, "Sobel", "thinning");
%blurred = conv2(edges, 1/9*ones(3, 3)) > 0;
img = bwmorph(edges,'thin', inf);

% ADJUST THE SECOND VALUE FOR DETAIL - larger value; less details
BW2 = bwareaopen(img, 10, 26); % cleans up the image by getting rid of small clusters

figure(2); imshow(BW2); title("BW Edge Image");

BW2 = padarray(BW2, [200 200]);
BW2 = fliplr(BW2');
[n, m]=size(BW2);
B=BW2; % Make a copy of the image (this will have parts erased as we go)
% Define 8 basic directions the Etch-a-sketch cursor can move
dirs=[0 1;0 -1;1 0;-1 0;1 1;-1 1;-1 -1;1 -1;2 1;-1 2;-2 -1;1 -2;1 2;-2 1;-1 -2;2 -1];
vec=[]; % This stores up the list of vectors which is the image in vector form
for i=1:n
    for j=1:m % Go through pixel by pixel to find a white one
        if B(i,j)==1
            for d=1:8 % From the white pixel, move in all 8 directions
                steps=1;
                while B(i+steps*dirs(d,1),j+steps*dirs(d,2))==1
                    steps=steps+1; % Count up how many steps you could go
                end                % while still being in a white space
                ss(d)=steps; % Store up the number of steps you could go
            end
            bestd=find(max(ss)==ss); % The direction where you could shoot out the farthest
            if max(ss)>=2 % Only keep vectors that are at least 4 long
                % Store the line segment as a start i j and an end i j 
                vec=[vec;[i j i+max(ss)*dirs(bestd(1),1) j+max(ss)*dirs(bestd(1),2)]];
                steps=1;
                for steps=0:max(ss) % Go back and erase the white pixels along the last stored line
                    B(i+steps*dirs(bestd(1),1),j+steps*dirs(bestd(1),2))=0;
                    r90=[0 1;-1 0]*dirs(bestd(1),:)'; % Erase 5 pixels to the right of the path
                    B(i+steps*dirs(bestd(1),1)+r90(1),j+steps*dirs(bestd(1),2)+r90(2))=0;
                    B(i+steps*dirs(bestd(1),1)+2*r90(1),j+steps*dirs(bestd(1),2)+2*r90(2))=0;
                    B(i+steps*dirs(bestd(1),1)+3*r90(1),j+steps*dirs(bestd(1),2)+3*r90(2))=0;
                    B(i+steps*dirs(bestd(1),1)+4*r90(1),j+steps*dirs(bestd(1),2)+4*r90(2))=0;
                    B(i+steps*dirs(bestd(1),1)+5*r90(1),j+steps*dirs(bestd(1),2)+5*r90(2))=0;
                    r90=[0 -1;1 0]*dirs(bestd(1),:)'; % Erase 5 pixels to the left of the path
                    B(i+steps*dirs(bestd(1),1)+r90(1),j+steps*dirs(bestd(1),2)+r90(2))=0;
                    B(i+steps*dirs(bestd(1),1)+2*r90(1),j+steps*dirs(bestd(1),2)+2*r90(2))=0;
                    B(i+steps*dirs(bestd(1),1)+3*r90(1),j+steps*dirs(bestd(1),2)+3*r90(2))=0;
                    B(i+steps*dirs(bestd(1),1)+4*r90(1),j+steps*dirs(bestd(1),2)+4*r90(2))=0;
                    B(i+steps*dirs(bestd(1),1)+5*r90(1),j+steps*dirs(bestd(1),2)+5*r90(2))=0;
                end
            end
            % Remember where your line segment ended
            lasti=i+max(ss)*dirs(bestd(1),1); 
            lastj=j+max(ss)*dirs(bestd(1),2);
            % Try making a chain of additional line segments
            for chain=1:10     
                for d=1:8    % Do the same thing you did before except starting from 
                    steps=1; % lasti,lastj instead of i,j
                    while B(lasti+steps*dirs(d,1),lastj+steps*dirs(d,2))==1
                        steps=steps+1;
                    end
                    ss(d)=steps;
                end
                bestd=find(max(ss)==ss);
                if max(ss)>=2
                    vec=[vec;[lasti lastj lasti+max(ss)*dirs(bestd(1),1) lastj+max(ss)*dirs(bestd(1),2)]];
                    steps=1;
                    for steps=0:max(ss) 
                        B(lasti+steps*dirs(bestd(1),1),lastj+steps*dirs(bestd(1),2))=0;
                        r90=[0 1;-1 0]*dirs(bestd(1),:)';
                        B(lasti+steps*dirs(bestd(1),1)+r90(1),lastj+steps*dirs(bestd(1),2)+r90(2))=0;
                        B(lasti+steps*dirs(bestd(1),1)+2*r90(1),lastj+steps*dirs(bestd(1),2)+2*r90(2))=0;
                        B(lasti+steps*dirs(bestd(1),1)+3*r90(1),lastj+steps*dirs(bestd(1),2)+3*r90(2))=0;
                        B(lasti+steps*dirs(bestd(1),1)+4*r90(1),lastj+steps*dirs(bestd(1),2)+4*r90(2))=0;
                        B(lasti+steps*dirs(bestd(1),1)+5*r90(1),lastj+steps*dirs(bestd(1),2)+5*r90(2))=0;
                        r90=[0 -1;1 0]*dirs(bestd(1),:)';
                        B(lasti+steps*dirs(bestd(1),1)+r90(1),lastj+steps*dirs(bestd(1),2)+r90(2))=0;
                        B(lasti+steps*dirs(bestd(1),1)+2*r90(1),lastj+steps*dirs(bestd(1),2)+2*r90(2))=0;
                        B(lasti+steps*dirs(bestd(1),1)+3*r90(1),lastj+steps*dirs(bestd(1),2)+3*r90(2))=0;
                        B(lasti+steps*dirs(bestd(1),1)+4*r90(1),lastj+steps*dirs(bestd(1),2)+4*r90(2))=0;
                        B(lasti+steps*dirs(bestd(1),1)+5*r90(1),lastj+steps*dirs(bestd(1),2)+5*r90(2))=0;
                    end
                end
                lasti=lasti+max(ss)*dirs(bestd(1),1); % Remember that last  
                lastj=lastj+max(ss)*dirs(bestd(1),2); % end of the chain
            end
        end
    end
end


figure(4);
axis equal
hold on
% expresses each vector as a set of two endpoints p1, p2
p1 = vec(:, 1:2); p2 = vec(:,3:4);
% sets the starting query point pq = p1(1,:)
pq = p1(1,:); n = 1; p1(n,:) = []; p2(n,:) = [];
h = animatedline;
title("Etch a Sketch")

while ~isempty(p1) && ~isempty(p2)
    
    % find the nearest point to the query point in each list
    nearest1 = dsearchn(p1(:,:), pq);
    nearest2 = dsearchn(p2(:,:), pq);
    
    % if the point is closest to a point in 1
    if norm(pq - p1(nearest1,:)) < norm(pq - p2(nearest2, :))
         n = nearest1;
         % plots a connecting point between pq and p1 (fills in the gap)
         addpoints(h, [pq(1) p1(n, 1)], [pq(2) p1(n, 2)]);
         % draws a line connecting p1 and p2 (the vector found in the
         % second part of the code)
         addpoints(h, [p1(n, 1) p2(n,1)], [p1(n,2) p2(n,2)]);
         % sets pq to the new endpoint and deletes already used points from
         % the matrices
         pq = p2(n,:); p1(n,:) = []; p2(n,:) = [];
         
    else % if the point is closest to a point in 2
         n = nearest2;
         addpoints(h, [pq(1) p2(n, 1)], [pq(2) p2(n, 2)]);
         addpoints(h, [p2(n,1) p1(n,1)], [p2(n,2) p1(n,2)]);
         pq = p1(n,:); p1(n,:) = []; p2(n,:) = []; 
    end
    % animates the etch a sketch. spicy
    % drawnow limitrate makes the drawing a lot faster
    drawnow  limitrate
end
