% Pareto front determination
% Weak: not worse in any criterium
% Minimization: Assumes minimum is better
%
% Example:
%   rawFits=rand(100,2); 
%   [bestFits, nonDomInds] = myParetoFront(rawFits,'plot');
%   rawFits=rand(100,3); 
%   [bestFits, nonDomInds] = myParetoFront(rawFits,'plot');
%
% Germán Augusto Ramírez, Nov 2019

function [bestFits, nonDomInds] = myParetoFront(rawFits,varargin)
    
    if ~isempty(varargin)
        if varargin{1} == 'plot'
            drawFront = true;
        end
    end
    ndims=size(rawFits,2);
    
    nonDomInds=[];
    for cont1 = 1:size(rawFits,1)
        dominated=false;
        for cont2 = 1:size(rawFits,1)
            if all(rawFits(cont1,:) > rawFits(cont2,:))
                dominated=true;
            end
        end
        if ~dominated
            nonDomInds = [nonDomInds; cont1];
        end
    end
    bestFits = rawFits(nonDomInds,:);

    if drawFront
        if ndims == 2 
            figure, scatter(rawFits(:,1),rawFits(:,2));
            [~,idx] = sort(bestFits(:,1));
            hold on; plot(bestFits(idx,1),bestFits(idx,2),'rx-'); 
            grid on, legend({'Data','Pareto Front'}); 
            legend({'Raw cost','Pareto front'}); 
            xlabel('Cost_x'); ylabel('Cost_y');
%             xx = linspace(min(bestFits(:,1)),max(bestFits(:,1)));
%             yy = interpn(sort(bestFits(idx,1)),bestFits(idx,2),xx,'spline'); % This can go realy bad
%             hold on, plot(xx,yy); 

        elseif ndims == 3
%             figure, plot3(rawFits(:,1),rawFits(:,2),rawFits(:,3),'bo');
%             hold on, plot3(bestFits(:,1),bestFits(:,2),bestFits(:,3),'r*');
%             xlabel('Cost_{SL}'); ylabel('Cost_{\theta}'); zlabel('Cost_{Icurr}');legend({'Data','Pareto Front'});

            F = scatteredInterpolant(bestFits(:,1),bestFits(:,2),bestFits(:,3),'linear','none');
            F.Method = 'linear';

            xx = linspace(min(bestFits(:,1)),max(bestFits(:,1)));
            yy = linspace(min(bestFits(:,2)),max(bestFits(:,2)));
            [XX,YY] = meshgrid(xx,yy);
            ZZ = F(XX,YY);

            figure
            surf(XX,YY,ZZ,'LineStyle','none'); % can also be mesh
            hold on, plot3(bestFits(:,1),bestFits(:,2),bestFits(:,3),'r*');
            hold on, plot3(rawFits(:,1),rawFits(:,2),rawFits(:,3),'bo');
            legend({'Raw cost','Pareto front'}); 
            xlabel('Cost_x'); ylabel('Cost_y'); zlabel('Cost_z')
        else
            disp('Im sorry, dont know how to plot that');
        end
    end
end