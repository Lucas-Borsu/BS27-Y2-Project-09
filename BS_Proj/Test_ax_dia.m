function [int_diam, ax_diam] = Test_ax_dia (int_dia, time_step)


    m= Cullen2018CortexAxonJPNlocalized_Kv12();
    m.intn.geo.diam.value.ref= int_dia;
    m.sim.dt.value=time_step;
    [membrane_pot, internode_length, time_vector] = ModelJPN_MTR(m);

    int_diam= m.intn.geo.diam.value.ref;
    ax_diam= m.node.geo.diam.value.ref;
    
end