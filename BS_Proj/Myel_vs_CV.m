
function [membrane_pot, internode_length, time_vector, myelin_thick] = Myel_vs_CV (n_myel, time_step)
    arguments
        n_myel = 20
        time_step = 1
    end

    n = n_myel;
    def_myel= 0.5894;
    idx75= round(n*0.75);
    max_val= def_myel*(n-1)/(idx75-1);

    myel_thick= linspace (0, max_val,n);

    node1 = 15;
    node2 = 35;

    membrane_pot = cell(1,n);
    internode_length = cell(1,n);
    time_vector = cell(1,n);
    myelin_thick = cell(1,n);
    cv= cell(1,n);

    for i = 1:n
        myelin_thick{i} = myel_thick(i);
      
        m= Cullen2018CortexAxonJPNlocalized_Kv12();
        m.sim.dt.value=time_step;

        m = UpdateInternodeSegmentDiameter(m, myelin_thick{i});
        
        [membrane_pot{i}, internode_length{i}, time_vector{i}] = ModelJPN_MTR(m);
        

        [~, idx1] = max(membrane_pot{i}(:,node1));          % ~ was max1
        [~, idx2] = max(membrane_pot{i}(:,node2));          % ~ was max2
        dist_um=sum(internode_length{i}(node1:node2,1));
        cv{i}= dist_um / ((idx2-idx1));

        %figure (1)
        %subplot(n/2,n/2,i)
        %plot(time_vector{i}, membrane_pot{i},'color',[0 i/n i/n])
        %title(['Voltage over time for myelin thickness = ' num2str(myelin_thick{i}) ' µm'])
        %xlabel('Time (ms)')
        %ylabel('Voltage (mV)')
    end

    x = cell2mat(myelin_thick);
    y = cell2mat(cv);
    
    figure(2)
    % 1. Create the line plot
    plot(x, y, '-o', 'LineWidth', 2, 'MarkerFaceColor', 'b', 'MarkerSize', 6)
    grid on
    
    xlabel('Myelin Thickness (µm)')
    ylabel('Conduction Velocity (µm/ms)')
    title('Conduction Velocity vs. Myelin Thickness')
    
    % 2. Add a text label precisely at the default myelin thickness (idx75)
    x_target = x(idx75);
    y_target = y(idx75);
    
    % text(x_pos, y_pos, 'String') - adjusted slightly so it doesn't overlap the data point
    text(x_target, y_target, ' \leftarrow base value = 0.5894', ...
        'FontSize', 10, 'FontWeight', 'bold', 'HorizontalAlignment', 'left');
end
