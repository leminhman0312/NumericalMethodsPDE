clear all
clc
imax = 3;
jmax = imax;
a = [1,2,3;4,5,6;7,8,9];


%shiftX
aplusX = zeros(imax,jmax);
% for i = 1:imax
%     for j = 1:jmax-1
%         aplusX(i,j) = a(i,j+1);
%     end
% end


% 
% aplusX(:,jmax) = a(:,1);
% 
% 
aplusX = circshift(a,-1,2);
aminusX = circshift(a,1,2);
% 
% 
% %shiftY 
% aplusY = circshift(a,1,1);
% 
% aminusY = zeros(imax,jmax);
% for j = 1:jmax
%       for i = 1:imax-1
%         aminusY(i,j) = a(i+1,j);
%       end
% end
%             
% aminusY(imax,:) = a(1,:);
% 
% 
% test1 = circshift(a,1,2);
% test2 = circshift(a,2,1);







