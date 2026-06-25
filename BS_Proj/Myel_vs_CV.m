
function [membrane_pot, internode_length, time_vector, myelin_thick] = Myel_vs_CV (n_myel, time_step)
    arguments
        n_myel = 25
        time_step = 1
    end

    n = n_myel;
    def_myel= 0.5894;
    idx60= round(n*0.60);

    step_size = def_myel/idx60;
    min_val= step_size;
    max_val = min_val + step_size * (n - 1);

    myel_thick= linspace (min_val, max_val,n);

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


        row_difference = idx2 - idx1;
        time_vec = time_vector{i}(idx2) - time_vector{i}(idx1);

        disp(['Row steps between peaks: ', num2str(row_difference)])
        disp(['Distance in time: ', num2str(time_vec)])


        cv{i} = dist_um / (time_vector{i}(idx2) - time_vector{i}(idx1));

        %figure (1)
        %subplot(n/2,n/2,i)
        %plot(time_vector{i}, membrane_pot{i},'color',[0 i/n i/n])
        %title(['Voltage over time for myelin thickness = ' num2str(myelin_thick{i}) ' µm'])
        %xlabel('Time (ms)')
        %ylabel('Voltage (mV)')
    end

    x = cell2mat(myelin_thick);
    y = cell2mat(cv);
    
    figure(1)
    % 1. Create the line plot
    plot(x, y, '-o', 'LineWidth', 2, 'MarkerFaceColor', 'b', 'MarkerSize', 6)
    grid on
    
    xlabel('Myelin Thickness (µm)')
    ylabel('Conduction Velocity (m/s)')
    title('Conduction Velocity vs. Myelin Thickness')
    
    ax = gca; % Get current axes
    
    % Get the standard numeric tick values MATLAB automatically picked
    default_ticks = ax.XTick;

    % 1. Create a clean list of ticks including your exact base value
    base_val = x(idx60);
    custom_ticks = unique(sort([default_ticks, base_val]));
    
    % 2. Find and remove any default tick that is too close to the base value
    % We use a small threshold (e.g., 0.03) to catch 0.6 without losing other numbers
    too_close = abs(custom_ticks - base_val) < 0.03 & (custom_ticks ~= base_val);
    custom_ticks(too_close) = []; % Deletes the overlapping tick (like 0.6)
    
    % 3. Apply the filtered ticks back to the axis
    ax.XTick = custom_ticks;
    
    % 4. Create and apply the string labels
    tick_labels = string(custom_ticks);
    base_tick_idx = find(custom_ticks == base_val);
    tick_labels(base_tick_idx) = "base value (" + string(base_val) + ")";
    ax.XTickLabel = tick_labels;
end
