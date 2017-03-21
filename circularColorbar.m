% circular colorbar

x=ones(1,360);
h=pie(x(1:2:end));
colormap hsv
for hc=1:2:length(h) % segments
    set(h(hc),'EdgeColor','none');
end
for hc=2:2:length(h) % texts
    if mod(hc,10)~=0
        delete(h(hc)); % delete some texts
    else
        set(h(hc),'string',num2str(hc),'fontSize',15); % set new texts
    end
end