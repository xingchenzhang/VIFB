function setfig(f)

% function setfig(f)
% Set figure dimensions for correct display of text to agree with
% print -deps.

figure(f)
set(gcf,'position',[220    40   790   526],'DefaultTextFontSize',14);
return

