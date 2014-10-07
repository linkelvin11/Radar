
tau= input('Enter pulse length (microsec): '); 
PRI= input('Enter pulse repetition interval (microsec): ');

CPI= input('Enter CPI microseconds): ');

tau=tau*1e-6;
PRI=PRI*1e-6;
CPI=CPI*1e-6;

%initialize hardcoded values
f=10e9;%x-band
c=3e8;

lambda=c/f;
range_res=1/CPI;
%range_res=tau*c/2
vel_res=lambda/2;
%vel_res=lambda*2
PRF=1/PRI;

%POSSIBLE CHANGING VALUES BASED ON CLUTTER
vel_width = round(20/vel_res);
%clutter blinds things within +/- 20 m/s 
range_width = round(((tau*c)+20)/range_res);
%clutter interferes with thing within 20 m


y_max=10*c*PRI/2;
x_max=4*PRF*lambda/2;

y=0:range_res:y_max;
x=0:vel_res:x_max;

Y=length(y);
X=length(x);

area=ones(Y,X);

area(:,1:vel_width)=0;


v1=abs((f/(f+PRF)-1)*c);
v1_index=round(v1/vel_res);
n=2;
while (v1_index<X)
    if((v1_index+vel_width)<X)
        area(:,(v1_index-vel_width):(v1_index+vel_width))=0;
    else
        area(:,(v1_index-vel_width):X)=0;
    end
    v1=abs((f/(f+(n*PRF))-1)*c);
    v1_index=round(v1/vel_res);
    n=n+1;
end

start = 1;
while (start < Y)
     if((start+range_width)<Y)
            area(start:(start+range_width),:)=0;
     else
         area(start:Y,:)=0;
     end
     start = start+round((PRI*c)/range_res);
end



area=area*255;

%y=y(length(y):-1:1);
y=-y;
figure
image(x,y,area);
colormap(gray(256))



