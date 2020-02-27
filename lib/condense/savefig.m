function savefig(obj)
    rez=1200; 
    f=gcf;
    figpos=getpixelposition(f);
    resolution=get(0,'ScreenPixelsPerInch'); 
    set(f,'paperunits','inches','papersize',figpos(3:4)/resolution,...
            'paperposition',[0 0 figpos(3:4)/resolution]);
    im = print('-RGBImage');
    [imind,cm] = rgb2ind(im,256);
    snapshot = char(strcat(obj.options.destination, 'step-', string(obj.iteration), '-clusters.png'));
    saveas(gcf, snapshot);
    filename = char(strcat(obj.options.destination, 'animated.gif'));

    if obj.iteration == 1
        imwrite(imind,cm,filename,'gif','Loopcount',inf,'DelayTime',0);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',1);
    end
    
    close all;
end