function [res] = plot_figures(T, THETA, PSI, HA, HP, SMA, ECC, k)

global R
 
if (k>1)
    
    SMA_plus = R/1000+375;
    
    SMA_minus = R/1000+300;

    i = 1;
    
    e = 0;
    
    while e <= 0.006
        
        hp_plus(i) = SMA_plus*(1-e)-R/1000;
        
        hp_minus(i) = SMA_minus*(1-e)-R/1000;
        
        ha_plus(i) = 2*(SMA_plus - R/1000)-hp_plus(i);
        
        ha_minus(i) = 2*(SMA_minus - R/1000)-hp_minus(i);
        
        e = e + 0.0001;
        
        i = i+1;
        
    end

    figure;

    subplot(1,2,1);
    
    plot(HP,HA,'o'),grid on, xlabel('HP (km)'),ylabel('HA (km)')
    
    hold on, plot([min(hp_minus),min(hp_minus)],[min(ha_minus), max(ha_plus)],'r','linewidth',2)
    
    plot([min(hp_minus),max(hp_plus)],[max(ha_plus), max(ha_plus)],'--b','linewidth',2)
    
    plot(hp_plus,ha_plus,'c','linewidth',2)
    
    plot(hp_minus,ha_minus,'c','linewidth',2)
    
    plot([300,375],[300,375],'k','linewidth',2)

    subplot(1,2,2);
    
    SMA_plus = SMA_plus - R/1000;
    
    SMA_minus = SMA_minus - R/1000;
    
    plot(ECC,SMA - R/1000,'o'),grid on, xlabel('Eccentricity'),ylabel('SMA Altitude (km)'), hold on
    
    plot([0,0.006],[SMA_plus , SMA_plus],'--r','linewidth',2)
    
    plot([0,0.006],[SMA_minus , SMA_minus],'--r','linewidth',2)
    
    plot([0,0],[SMA_minus, SMA_plus],'--k','linewidth',2)
    
    plot([0.006,0.006],[SMA_minus, SMA_plus],'--k','linewidth',2)
    
end


if (k == 1)
    
    figure, plot(T,THETA*180/pi,'linewidth',2), grid on
    
    ylabel('Theta (deg)') 
    
    figure, plot(T,PSI*180/pi,'linewidth',2), grid on
    
    ylabel('Psi (deg)')
    
end

res = 0;

end