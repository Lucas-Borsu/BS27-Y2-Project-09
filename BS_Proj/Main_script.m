
function [membrane_pot, internode_length, time_vector, myelin_thick] = Main_script (myelin_thickness,g_ratio, time_step)
    n = length(myelin_thickness);
        
    if n ~= length(g_ratio)
        error('Input needs to be of the same length')
    end
    
    node1 = 15;
    node2 = 35;

    membrane_pot = cell(1,n);
    internode_length = cell(1,n);
    time_vector = cell(1,n);
    myelin_thick = cell(1,n);
    gratio = cell(1,n);
    cv= cell(1,n);

    for i = 1:n
        myelin_thick{i} = myelin_thickness{i};
        gratio{i} = g_ratio{i};

        m= Cullen2018CortexAxonJPNlocalized_Kv12();
        m.sim.dt.value=time_step;

        m = UpdateInternodeSegmentDiameter(m, myelin_thick{i});
        m = UpdateInternodeGRatio(m, gratio{i});
        
        [membrane_pot{i}, internode_length{i}, time_vector{i}] = ModelJPN_MTR(m);
        

        [max1, idx1] = max(membrane_pot{i}(:,node1));
        [max2, idx2] = max(membrane_pot{i}(:,node2));
        dist_um=sum(internode_length{i}(node1:node2,1));
        cv{i}= dist_um / ((idx2-idx1));

        figure (1)
        subplot(3,3,i)
        plot(time_vector{i}, membrane_pot{i},'color',[0 i/n i/n])
        title(['Voltage over time for myelin thickness = ' num2str(myelin_thickness{i}) ' µm'])
        xlabel('Time (ms)')
        ylabel('Voltage (mV)')
    end

    x = cell2mat(myelin_thick);
    y = cell2mat(cv);
    
    figure(2)
    b = bar(categorical(string(x)), y);
    
    %ylim([min(y)*0.98 max(y)*1.02]);

    b.FaceColor = 'flat';
    b.CData = parula(n);
    
    xlabel('Myelin Thickness (µm)')
    ylabel('Conduction Velocity (µm/ms)')
end
