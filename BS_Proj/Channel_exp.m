function Channel_exp (myelin_thickness, time_step)
    
    node1 = 15;
    node2 = 35;

    membrane_pot = cell(1,2);
    internode_length = cell(1, 2);
    time_vector = cell(1, 2);
    cv = cell(1, 2);

    Kv12_exp= Cullen2018CortexAxonJPNlocalized_Kv12_exp();
    Kv12_exp.sim.dt.value=time_step;
    
    Kv12_exp = UpdateInternodeGRatio(Kv12_exp, 0.999, [], [1, 2, 3, 4]);
            
    Kv12_exp.channels(4).location = [Kv12_exp.geo.juxtaparanodeSegments{:}]';
    Kv12_exp.active.segments = unique(cat(1, Kv12_exp.channels(:).location));

    [membrane_pot{1}, internode_length{1}, time_vector{1}] = ModelJPN_MTR(Kv12_exp);


    m= Cullen2018CortexAxonJPNlocalized_Kv12();
    m.sim.dt.value=time_step;

    m = UpdateInternodeSegmentDiameter(m, myelin_thickness);
    
        
    [membrane_pot{2}, internode_length{2}, time_vector{2}] = ModelJPN_MTR(m);

    [~, idx1] = max(membrane_pot{1}(:,node1));          % ~ was max1
    [~, idx2] = max(membrane_pot{1}(:,node2));          % ~ was max2
    dist_um=sum(internode_length{1}(node1:node2,1));
    cv{1}= dist_um / ((idx2-idx1));

    [~, idx1] = max(membrane_pot{2}(:,node1));          % ~ was max1
    [~, idx2] = max(membrane_pot{2}(:,node2));          % ~ was max2
    dist_um=sum(internode_length{2}(node1:node2,1));
    cv{2}= dist_um / ((idx2-idx1));

    

    cv_data = cell2mat(cv);
    
    % 2. Open a new figure window
    figure('Name', 'Conduction Velocity Comparison', 'NumberTitle', 'off');
    
    % 3. Define categorical names for the X-axis bars
    categories = categorical({'Stripped Myelin (g=1.0)', 'Modified Thickness'});
    % Reorder categories so MATLAB doesn't sort them alphabetically
    categories = reordercats(categories, {'Stripped Myelin (g=1.0)', 'Modified Thickness'});
    
    % 4. Create the bar chart
    b = bar(categories, cv_data);
    
    % 5. Aesthetic customizations (Colors, Labels, Grid)
    b.FaceColor = 'flat';
    b.CData(1,:) = [0.8500 0.3250 0.0980]; % Red/Orange for stripped myelin
    b.CData(2,:) = [0 0.4470 0.7410];      % Blue for modified thickness
    
    ylabel('Conduction Velocity (µm/sample index)');
    title('Conduction Velocity Comparison between Myelin Conditions');
    grid on;
    
    % Set y-axis limits dynamically with a 5% margin above the max value
    ylim([0, max(cv_data) * 1.05]);
    
    figure(2)
    plot(time_vector{2}', membrane_pot{2}(:, 35), 'b', 'LineWidth', 2.5);  hold on;
    plot(time_vector{1}', membrane_pot{1}(:, 35), 'r');
    
    
end