
function [Source,Receiver] = Point_generation(Lx,Ly,Lz,flag)
    % This
    % Input:
    %  -num:    Source/Receiver amount 
    %  -Lx:     Length of the room  
    %  -Ly:     Width of the room
    %  -Lz:     Height of the room
    %  -flag:   Whether show the figure of source/receiver distribution 
    %  
    %
    % Output:
    %  -rho:    Density of air, kg/m3
    %  -c0:     Zero-frequency speed of sound in air, m/s
   
%% Source grid
nx = 3;
ny = 1;
nz = 2;

x1 = zeros(nx*ny*nz,1);
y1 = x1;
z1 = x1;

% ny = 1
for ix = 0:nx-1
    for iz = 1:nz
        x1(ix*nz+iz) = (rand(1)+ix)*Lx/nx;
        y1(ix*nz+iz) =  rand(1)*Ly;
        z1(ix*nz+iz) = ((iz-1)+rand(1))*Lz/nz;
    end
end


%% Receiver grid

nx = 6;
ny = 6;
nz = 6;

x2 = zeros(nx*ny*nz,1);
y2 = x2;
z2 = x2;


for ix = 0:nx-1
    for iy = 0:ny-1
        for iz = 1:nz
            x2(iy*(nx*nz)+ix*nz+iz) = (rand(1)+ix)*Lx/nx;
            y2(iy*(nx*nz)+ix*nz+iz) =  (iy+rand(1))*Ly/ny;
            z2(iy*(nx*nz)+ix*nz+iz) = ((iz-1)+rand(1))*Lz/nz;
        end
    end
end

Source = [x1,y1,z1];
Receiver = [x2,y2,z2];

    
    %% Source point ranking

    % According to the numbering rule in COMSOL
    % 1. First compare x coordinate, The smaller the x, the smaller the number
    % 2. If x1 = x2, then compare y, the smaller the y, the smaller the number
    % 3. If x1 = x2, y1 = y2, then compare z

    R  = sortrows(Source,1);

    countx = hist(R(:,1),unique(R(:,1))); % See if each x is the unique one
    
    k = find(countx-1); % Find the none 1 elemnt (none unique x)
    
    if k ~=0
        for i = 1:numel(k) % if some x is non unique, then the order these point order is determined by y
            Startx = sum(countx(1:k(i)-1))+1;
            Endx = sum(countx(1:k(i)-1))+countx(k(i));
            tempx = R(Startx:Endx,:);
            tempx = sortrows(tempx,2);
                county = hist(tempx (:,2),unique(tempx (:,2))); % See if each y is the unique one
                h = find(county-1); % Find the none 1 elemnt (none unique y)
                if h ~=0
                    for j = 1:numel(h) % if some x is non unique, then the order these point order is determined by y
                        Starty = sum(county(1:h(j)-1))+1;
                        Endy = sum(county(1:h(j)-1))+county(h(j));
                        tempy = tempx(Starty:Endy,:);
                        tempy = sortrows(tempy,3);
                        tempx(Starty:Endy,:) = tempy;
                    end
                end
            R(Startx:Endx,:) = tempx;
        
        end
    end
    
    Source = R;
    clear R

    %% Receiver point ranking

    % According to the numbering rule in COMSOL
    % 1. First compare x coordinate, The smaller the x, the smaller the number
    % 2. If x1 = x2, then compare y, the smaller the y, the smaller the number
    % 3. If x1 = x2, y1 = y2, then compare z

    R  = sortrows(Receiver,1);

    countx = hist(R(:,1),unique(R(:,1))); % See if each x is the unique one
    
    k = find(countx-1); % Find the none 1 elemnt (none unique x)
    
    if k ~=0
        for i = 1:numel(k) % if some x is non unique, then the order these point order is determined by y
            Startx = sum(countx(1:k(i)-1))+1;
            Endx = sum(countx(1:k(i)-1))+countx(k(i));
            tempx = R(Startx:Endx,:);
            tempx = sortrows(tempx,2);
                county = hist(tempx (:,2),unique(tempx (:,2))); % See if each y is the unique one
                h = find(county-1); % Find the none 1 elemnt (none unique y)
                if h ~=0
                    for j = 1:numel(h) % if some x is non unique, then the order these point order is determined by y
                        Starty = sum(county(1:h(j)-1))+1;
                        Endy = sum(county(1:h(j)-1))+county(h(j));
                        tempy = tempx(Starty:Endy,:);
                        tempy = sortrows(tempy,3);
                        tempx(Starty:Endy,:) = tempy;
                    end
                end
            R(Startx:Endx,:) = tempx;
        
        end
    end
    Receiver = R;



     %% Cube
     if nargin == 4 && flag == 1
        % show the points distribution

        x_inner=Lx/2;     %length
        y_inner=Ly;     %width
        z_inner=Lz;     %height
        % Define the intersection position of each plane and axis ;
        x_outer=x_inner+Lx;
        
        x_inner1=0;     
        x_inner2=x_inner;  
        x_outer1=x_outer/2+x_inner/2;     
        x_outer2=0;
        
        
        y_inner1=y_inner;     
        y_inner2=0;
        y_outer1=y_inner;     
        y_outer2=0; 
        
        
        z_inner1=z_inner;     
        z_inner2=0;     
        z_outer1=z_inner;     
        z_outer2=0;
        figure
        % Draw the inner cube box （ Inspection quality block ）;
        % Vertex matrix ;
        vertex_matrix=[x_inner2 y_inner2 z_inner2;x_inner1 y_inner2 z_inner2;x_inner1 y_inner1 z_inner2;x_inner2 y_inner1 z_inner2;
                       x_inner2 y_inner2 z_inner1;x_inner1 y_inner2 z_inner1;x_inner1 y_inner1 z_inner1;x_inner2 y_inner1 z_inner1];
        
        % Connection matrix ： The value in each row of the connection matrix represents the row mark of the vertex matrix ;
        face_matrix=[1 2 6 5;2 3 7 6;3 4 8 7;
                     4 1 5 8;1 2 3 4;5 6 7 8];
        patch('Vertices',vertex_matrix,'Faces',face_matrix,'FaceVertexCData',hsv(8),'FaceColor','none')
        
        
        % Draw the outer cube frame （ Electrode cage ）;
        % Vertex matrix ;
        vertex_matrix=[x_outer2 y_outer2 z_outer2;x_outer1 y_outer2 z_outer2;x_outer1 y_outer1 z_outer2;x_outer2 y_outer1 z_outer2;
                       x_outer2 y_outer2 z_outer1;x_outer1 y_outer2 z_outer1;x_outer1 y_outer1 z_outer1;x_outer2 y_outer1 z_outer1];
        
        % Connection matrix ： The value in each row of the connection matrix represents the row mark of the vertex matrix ;
        face_matrix=[1 2 6 5;2 3 7 6;3 4 8 7;
                     4 1 5 8;1 2 3 4;5 6 7 8];
        patch('Vertices',vertex_matrix,'Faces',face_matrix,'FaceVertexCData',hsv(8),'FaceColor','none')
        hold on
        
        view(3);
        scatter3(x1, y1, z1, 'MarkerFaceColor' , 'b' , 'MarkerEdgeColor' , 'b' );
        scatter3(x2, y2, z2, 'MarkerFaceColor' , 'r' , 'MarkerEdgeColor' , 'r' );
     end 
end
