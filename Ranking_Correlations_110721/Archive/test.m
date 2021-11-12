figure
for num = 1:4
    if num == 1
        clr = [0, 0.4470, 0.7410];
    elseif num == 2
        clr = [.93,.69,.13];
    elseif num == 3
        clr = [.85,.33,.1];
    elseif num == 4
        clr = [.64,.08,.18];
    end
    
    scatplot=scatter(0,0,100,'o','MarkerFaceColor',clr,'MarkerEdgeColor',clr); 
    % Set property MarkerFaceAlpha and MarkerEdgeAlpha to <1.0
    scatplot.MarkerFaceAlpha = .3;
    hold on
end

legend('TD','ASD Severity 1','ASD Severity 2','ASD Severity 3','Orientation','horizontal')