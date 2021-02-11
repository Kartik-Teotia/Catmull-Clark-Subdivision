close all
clear;
[V,F] = loadmesh('mushroom.off');
V = V';
F = F';


figure(1)
plotquadmesh(V,F);


for it = 1:3
input_faces = F;
input_points = V;
[V,F]=CCS(input_faces,input_points);
figure(it+1)
plotquadmesh(V,F);

end

function [VV,FF] = CCS(input_faces,input_points)
input_faces = input_faces -1;
a=[0;0;1];
b=[0,0,0];
cs=sum_point(a,b);
face_points = get_face_point(input_points, (input_faces+1)');
edges_faces = get_edges_faces(input_points, (input_faces+1));
edge_points =  get_edge_points(input_points,edges_faces,face_points);
avg_face_points = get_avg_face_points(input_points, input_faces+1, face_points);
avg_mid_edges = get_avg_mid_edges(input_points,edges_faces);
points_faces = get_points_faces(input_points, (input_faces+1));
new_points = get_new_points(input_points,points_faces,avg_face_points',avg_mid_edges');
face_point_nums = [];
next_pointnum = length(new_points);
for face_point = face_points'
    
    new_points = [new_points ; face_point'];
    face_point_nums=[face_point_nums;next_pointnum];
    next_pointnum = next_pointnum + 1;
end

    edge_point_nums = [];
    edge_points=edge_points';
    for edgenum = 1:length(edges_faces)
        pointnum_1 = edges_faces(1,edgenum);
        pointnum_2 = edges_faces(2,edgenum);
        edge_point = edge_points(edgenum,:);
        new_points = [new_points; edge_point];
        edge_point_nums(pointnum_1, pointnum_2) = next_pointnum;
        next_pointnum = next_pointnum + 1;
    end
input_faces = input_faces + 1;    
new_faces = [];
for oldfacenum = 1: length(input_faces)
    oldface = input_faces(oldfacenum,:);
    if length(oldface) == 3
        a = oldface(1);
        b = oldface(2);
        c = oldface(3);
        face_point_abcd = face_point_nums(oldfacenum);
          w=switch_nums([a,b]);
    a1=w(1);
    b1=w(2);
    edge_point_ab  = edge_point_nums(a1,b1);
    x=switch_nums([c,a]);
    a2 = x(1);
    b2 = x(2);
    y=switch_nums([b,c]);
    a3 = y(1);
    b3 = y(2);
    edge_point_ca  = edge_point_nums(a2,b2);
    edge_point_bc  = edge_point_nums(a3,b3);
    %edge_point_cd  = edge_point_nums(a4,b4);
    new_faces = [new_faces; a,edge_point_ab,face_point_abcd,edge_point_ca];
    new_faces = [new_faces; b,edge_point_bc,face_point_abcd,edge_point_ab];
    new_faces = [new_faces; c,edge_point_ca,face_point_abcd,edge_point_bc];
    %new_faces = [new_faces; d,edge_point_da,face_point_abcd,edge_point_cd]
    else   
    a = oldface(1);
    b = oldface(2);
    c = oldface(3); 
    d = oldface(4); 
    face_point_abcd = face_point_nums(oldfacenum);
    w=switch_nums([a,b]);
    a1=w(1);
    b1=w(2);
    edge_point_ab  = edge_point_nums(a1,b1);
    x=switch_nums([d,a]);
    a2 = x(1);
    b2 = x(2);
    y=switch_nums([b,c]);
    a3 = y(1);
    b3 = y(2);
    z=switch_nums([c,d]);
    a4 = z(1);
    b4 = z(2);
    edge_point_da  = edge_point_nums(a2,b2);
    edge_point_bc  = edge_point_nums(a3,b3);
    edge_point_cd  = edge_point_nums(a4,b4);
    new_faces = [new_faces; a,edge_point_ab,face_point_abcd,edge_point_da];
    new_faces = [new_faces; b,edge_point_bc,face_point_abcd,edge_point_ab];
    new_faces = [new_faces; c,edge_point_cd,face_point_abcd,edge_point_bc];
    new_faces = [new_faces; d,edge_point_da,face_point_abcd,edge_point_cd];
    end
end
new_faces(:,1) = new_faces(:,1) -1;
new_faces = new_faces + 1;
VV = new_points;
FF = new_faces;
end
    
%edge_point_nums=sparse(edge_point_nums)  
function face_points  = get_face_point(input_points, input_faces)
num_dimensions = 3;
face_points = [];
for curr_face=input_faces
    
    face_point = [0,0,0];
    for curr_point_index = curr_face;


        curr_point = input_points(curr_point_index,:);
        face_point=sum(curr_point);
    
        num_points = length(curr_face);
        face_points=[face_points;face_point/length(curr_face)];
        
    end
    

end
end
function edges_centers = get_edges_faces(input_points, input_faces);
edges = [];
for facenum = 1 : length(input_faces)
    face = input_faces(facenum,:);
    num_points = length(face);
    for pointindex = 1:num_points
      
        if pointindex< num_points
            pointnum_1 = face(pointindex);
            pointnum_2 = face(pointindex+1);
        else
   
        pointnum_1 = face(pointindex);
        pointnum_2 = face(1);
        end
        if pointnum_1>pointnum_2
            temp  = pointnum_1;
            pointnum_1=pointnum_2;
            pointnum_2 = temp;
        end
     edges = [edges; [pointnum_1,pointnum_2,facenum]];
    end
end
edges = sortrows(edges);
num_edges = length(edges);
            eindex = 1;
            merged_edges = [];
            while eindex <= num_edges
                e1 = edges(eindex,:);
                if eindex < num_edges             
                e2 = edges(eindex + 1,:);
                
                if e1(1) == e2(1) && e1(2) == e2(2)
                    merged_edges = [merged_edges; [e1(1),e1(2),e1(3),e2(3)]];
                        eindex = eindex + 2;
                   
                 else
                merged_edges = [merged_edges; e1(1),e1(2),e1(3), NaN];
                eindex = eindex + 1;
                end
                else
                merged_edges = [merged_edges; e1(1), e1(2), e1(3), NaN];
                eindex = eindex + 1;
                end
                
            end
            edges_centers = [];
            for me = merged_edges'
                p1 = input_points(me(1),:);
                p2 = input_points(me(2),:);
                cp = center_point(p1,p2);
                edges_centers = [edges_centers, cat(1,me,[cp])];
            end 
end
function edge_points = get_edge_points(input_points, edges_faces, face_points)
edge_points = [];

for edge = edges_faces

        cp = edge(5:end);
    fp1 = face_points(edge(3),:);
    if isnan(edge(4))
        fp1 = fp2;
    else
       fp2 = face_points (edge(4),:);
    end
    cfp = center_point(fp1,fp2);
    edge_points = [edge_points, center_point(cp,cfp)];
   
end
end
function avg_face_points = get_avg_face_points(input_points, input_faces, face_points)
num_points = length(input_points);
temp_points = []; 
for pointnum = 1:num_points
    temp_points = [temp_points; 0,0,0,0];
end

for facenum = 1:length(input_faces)
    
    fp = face_points(facenum,:);
    for pointnum = input_faces(facenum,:)
        
        %tp = temp_points(pointnum,:)
        tp = temp_points(pointnum,1:3);
        
        temp_points(pointnum,1:3) = fp+tp;
        temp_points(pointnum,4) = temp_points(pointnum,4) + 1;
        
    end
 
end

    
   avg_face_points = [];
   for tp = temp_points'
       
       afp = tp(1:3)./ tp(4);
       avg_face_points = [avg_face_points, afp];
   end
   
end
function avg_mid_edges = get_avg_mid_edges(input_points,edges_faces)
num_points = length(input_points);
temp_points = [];
for pointnum = 1:num_points
    temp_points = [temp_points; 0,0,0,0];
end

for edge = [edges_faces];
    cp=edge(5:7);
    for pointnum = [edge(1),edge(2)]
    
        tp = temp_points(pointnum,1:3);
        temp_points(pointnum,1:3) = cp' + tp;
        temp_points(pointnum,4) = temp_points(pointnum,4) + 1;
    end
end
avg_mid_edges = [];
for tp = temp_points'
    
    ame = tp(1:3)./tp(4);
    avg_mid_edges = [avg_mid_edges, ame];
end
end
function points_faces = get_points_faces(input_points, input_faces)
num_points = length(input_points)

points_faces = [];
for pointnum = 1:num_points
    points_faces = [points_faces , 0];
end

for facenum = 1:length(input_faces)
    for pointnum = input_faces(facenum,:)
    points_faces(pointnum) = points_faces(pointnum) + 1;
    end
end
end
function new_points = get_new_points(input_points,points_faces,avg_face_points,avg_mid_edges)
new_points = [];
for pointnum = 1:length(input_points)
    n = points_faces(pointnum);
    m1 = (n-3)/n;
    m2= 1/n;
    m3 = 2/n;
    old_coords = input_points(pointnum,:);
    p1 = old_coords.*m1;
    afp=avg_face_points(pointnum,:);
    p2=afp.*m2;
    ame=avg_mid_edges(pointnum,:);
    p3=ame.*m3;
    p4=p1+p2;
    new_coords=p4+p3;
    new_points = [new_points; new_coords];
end
end
function point_nums = switch_nums(point_nums)
if point_nums(1)<point_nums(2);
    point_nums = point_nums;
else
    point_nums=[point_nums(2),point_nums(1)];
end
        
end
function sp =div_point(p,d)
sp=[]
for it=1:3
    sp = [sp, p(it)/d];
end
end
function sp= sum_point(p1,p2)
sp=[];
for it=1:3
    sp = [sp, p1(it) + p1(it)];
end
end    

function cp = center_point(p1,p2)
cp=[];
for it=1:3
    cp = [cp;(p1(it)+p2(it))/2];
end
end         

     