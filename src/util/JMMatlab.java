package util;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.ListIterator;

/*
 * tuas variáveis de entrada são PS, e o array TC_i
 */
public class JMMatlab {

    public static void main(String[] args) {

        List<Double> TC_i = new ArrayList<>();
        TC_i.addAll(Collections.nCopies(3, 1.d));


        JMMatlab jmMatlab = new JMMatlab();
        double lat = jmMatlab.latencia(300, TC_i);
        System.out.println(lat);
    }

    public double latencia(float PS, List<Double> TC_i) {

        double lats = 0.d;

        int h = 2; // % # hops;
        double TC_l = 1; // % cycles;
        double TS = 6; // %cycles;
        double TS_l = 3; // % cycles;

        // TC_i = [TC_l 1/(1-FR) ones(1,h-1) TC_l];
        List<Double> TS_i = new ArrayList<Double>();
        TS_i.add(PS - TC_l);
        TS_i.addAll(Collections.nCopies(TC_i.size() - 2, TS));
        TS_i.add((TS_l));

        // phi_out = 1./TC_i; % flits/cycle
        List<Double> phi_out = new ArrayList<Double>();
        for (final ListIterator<Double> i = TC_i.listIterator(); i.hasNext();) {
            final double element = i.next();
            phi_out.add(1 / element);
        }

        // phi_in = [1 phi_out(1:end-1)];
        List<Double> phi_in = new ArrayList<Double>();
        List<Double> phi_outCOPY = new ArrayList<Double>(phi_out);
        phi_outCOPY.remove(phi_out.size() - 1);
        phi_in.add(1.d);
        phi_in.addAll(phi_outCOPY);

        // V0_i = (TS_i + TC_i).*phi_in; % flits
        List<Double> V0_i = new ArrayList<>();
        for (int i = 0; i < TS_i.size(); i++) {
            V0_i.add((TS_i.get(i) + TC_i.get(i)) / phi_in.get(i));
        }

        // phi_out_eff = phi_out;
        List<Double> phi_out_eff = new ArrayList<>(phi_out);
        // for i = 2:length(phi_out_eff)
        for (int i = 1; i < phi_out_eff.size(); i++) {
            // if phi_out_eff(i) > phi_out_eff(i-1),
            if (phi_out_eff.get(i) > phi_out_eff.get(i - 1)) {
                // V0_i(i) = (TS_i(i) + TC_i(i))*phi_out_eff(i-1);
                V0_i.set(i, (TS_i.get(i) + TC_i.get(i)) * phi_out_eff.get(i - 1));

                // TE = V0_i(i)/(phi_out(i) - phi_out_eff(i-1));
                double TE = V0_i.get(i) / (phi_out.get(i) - phi_out_eff.get(i - 1));

                // T2 = (PS - V0_i(i) - TE*phi_out_eff(i-1))/(phi_out_eff(i-1));
                double T2 = (PS - V0_i.get(i) - TE * phi_out_eff.get(i - 1)) / (phi_out_eff.get(i - 1));

                // phi_out_eff(i) = (phi_out(i)*TE + phi_out_eff(i-1)*T2)/(TE +
                // T2);
                phi_out_eff.set(i, (phi_out.get(i) * TE + phi_out_eff.get(i - 1) * T2) / (TE + T2));
            }
        }

        // phi_in_eff = [1 phi_out_eff(1:end-1)];
        List<Double> phi_in_eff = new ArrayList<Double>();
        List<Double> phi_out_effCOPY = new ArrayList<Double>(phi_out_eff);
        phi_out_effCOPY.remove(phi_out_eff.size() - 1);
        phi_in_eff.add(1.d);
        phi_in_eff.addAll(phi_out_effCOPY);


        // tau_i = (PS - V0_i)./phi_in; //% cycles
        List<Double> tau_i = new ArrayList<>();
        for (int i = 0; i < V0_i.size(); i++) {
            tau_i.add((PS - V0_i.get(i)) / phi_in_eff.get(i));
        }

        //Phi_i = phi_in_eff - phi_out_eff;
        List<Double> Phi_i = new ArrayList<Double>();
        for (int i = 0; i < TS_i.size(); i++) {
            Phi_i.add(phi_in_eff.get(i) - phi_out_eff.get(i));
        }

        // V_tau_i = max([V0_i + Phi_i.*tau_i; zeros(size(V0_i))]);
        List<Double> V_tau_i = new ArrayList<Double>();
        for (int i = 0; i < Phi_i.size(); i++) {
            double operacao = V0_i.get(i) + (Phi_i.get(i) * tau_i.get(i));
            V_tau_i.add(operacao > 0 ? operacao : 0.d);
        }

        // V_tau_i(find(V_tau_i == 0)) = 1;
        Collections.replaceAll(V_tau_i, 0.d, 1.d);


        //   DT_i = V_tau_i./phi_out_eff;
        List<Double> DT_i = new ArrayList<Double>();
        for (int i = 0; i < V_tau_i.size(); i++) {
            DT_i.add(V_tau_i.get(i) / phi_out_eff.get(i));
        }

        // lats = [lats sum(DT_i)];
        for (Double double1 : DT_i) {
            lats += double1;
        }

        return lats;
    }

}