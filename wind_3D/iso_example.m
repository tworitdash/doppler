 func = @(x,y,z)(sin(2*x.^2+2*y.^2).*cos(2*z.^2));
 [x,y,z] = meshgrid(-1:0.1:1, -1:0.1:1, -1:0.1:1);
 v = func(x,y,z);
 figure
 p1 = patch(isosurface(x,y,z,v,0.5));
 hold on
 p2 = patch(isosurface(x,y,z,v,0.8));
 isonormals(x,y,z,v,p1);
 isonormals(x,y,z,v,p2);
 p1.FaceColor = 'red';
 p2.FaceColor = 'green';
 p1.EdgeColor = 'none';
 p2.EdgeColor = 'none';
 daspect([1,1,1])
 view(3); axis tight
 camlight 
 lighting gouraud