using Graphs, LinearAlgebra, Statistics, Plots, DelimitedFiles, StatsBase, SparseArrays

# --- 1. Entropia Estrutural (Ginestra Bianconi) ---
function bianconi_entropy(g)
    degrees = degree(g)
    N = nv(g)
    k_avg = mean(degrees)
    S = 0.0
    for i in 1:N
        for j in (i+1):N
            p_ij = (degrees[i] * degrees[j]) / (N * k_avg)
            p_ij = clamp(p_ij, 1e-10, 1 - 1e-10)
            S -= (p_ij * log(p_ij) + (1 - p_ij) * log(1 - p_ij))
        end
    end
    return S / N
end

# --- 2. Teoria da Informação (Shannon H e Mutual Info MI) ---
function calculate_info_metrics(phases)
    N, T = size(phases)
    h_bins = fit(Histogram, vec(phases[:, end]), -π:0.2:π).weights
    p_dist = h_bins ./ sum(h_bins)
    p_dist = p_dist[p_dist .> 0]
    shannon_h = -sum(p_dist .* log2.(p_dist))
    
    mi_total, count = 0.0, 0
    # Amostra restrita aos primeiros 20 nós para otimização do cálculo
    for i in 1:min(N, 20), j in (i+1):min(N, 20)
        h_i = fit(Histogram, phases[i, :], -π:0.2:π).weights .+ 1e-10
        h_j = fit(Histogram, phases[j, :], -π:0.2:π).weights .+ 1e-10
        h_ij = fit(Histogram, (phases[i, :], phases[j, :]), (-π:0.2:π, -π:0.2:π)).weights .+ 1e-10
        p_i, p_j, p_ij = h_i./sum(h_i), h_j./sum(h_j), h_ij./sum(h_ij)
        mi_total += sum(p_i .* log2.(p_i)) + sum(p_j .* log2.(p_j)) - sum(p_ij .* log2.(p_ij))
        count += 1
    end
    return shannon_h, abs(mi_total / count)
end

# --- 3. Solucionador Piezo-Duplex Fracionário ---
function solve_piezo_duplex_fractional(N, L1, L2, alpha, p_phys, t_steps, dt)
    η, κ_hys, Abw, β, γ, n_bw, f0, ω, σ_mech, σ_elec, σ_inter = p_phys
    cp = ones(t_steps)
    for j in 1:t_steps-1; cp[j+1] = cp[j] * (1 - (alpha + 1)/j); end

    x = zeros(N, t_steps); v = zeros(N, t_steps); z = zeros(N, t_steps)
    V_el = zeros(N, t_steps)
    
    x[:, 1] .= randn(N) * 0.05; V_el[:, 1] .= randn(N) * 0.05
    dt_a = dt^alpha
    
    # Parâmetros de modulação de fase (Broadband excitation)
    a0, b0 = 0.5, 0.5 

    for t in 1:t_steps-1
        m_x = zeros(N); m_v = zeros(N); m_z = zeros(N); m_V = zeros(N)
        if t > 1
            for j in 1:t-1
                c = cp[j+1]
                @. m_x += c * x[:, t-j]
                @. m_v += c * v[:, t-j]
                @. m_z += c * z[:, t-j]
                @. m_V += c * V_el[:, t-j]
            end
        end

        c_mech = -σ_mech * (L1 * x[:, t])
        c_elec = -σ_elec * (L2 * V_el[:, t])
        
        # Força de excitação com modulação de fase
        time_t = t * dt
        F_ext = f0 * cos(ω * time_t + a0 * sin(b0 * ω * time_t))

        # Atualização com travas físicas de saturação
        @. x[:, t+1] = clamp(v[:, t] * dt_a - m_x, -50.0, 50.0)
        
        v_new = @. (-2η*v[:, t] + 0.5*x[:, t]*(1 - x[:, t]^2) - κ_hys*z[:, t] + F_ext + c_mech + σ_inter*V_el[:, t]) * dt_a - m_v
        @. v[:, t+1] = clamp(v_new, -50.0, 50.0)
        
        z_new = @. (Abw*v[:, t] - β*abs(v[:, t])*z[:, t]*abs(z[:, t])^(n_bw-1) - γ*v[:, t]*abs(z[:, t])^n_bw) * dt_a - m_z
        @. z[:, t+1] = clamp(z_new, -10.0, 10.0)
        
        V_new = @. (-V_el[:, t] - σ_inter*v[:, t] + c_elec) * dt_a - m_V
        @. V_el[:, t+1] = clamp(V_new, -100.0, 100.0)
    end
    return x, v, V_el
end

# --- 4. Simulação em Lote (Batch Simulation) ---

N_values = [100, 200, 500, 1000, 1500, 2000]
alphas = [0.6, 0.8, 1.0]
sigma_inter_range = range(0.0, 1.0, length=100) # Alta resolução com 100 pontos
dt, t_steps = 0.001, 10000 # Tempo total aumentado para 10s e passo super-preciso

for N in N_values
    println("\n======================================")
    println("Iniciando bateria para N = $N")
    println("======================================")
    
    # Cria a pasta exclusiva para este tamanho de rede
    folder = "Resultados_N_$N"
    mkpath(folder)

    # Construção das Redes
    g_mech = watts_strogatz(N, 4, 0.1) 
    g_elec = erdos_renyi(N, 4/N)      
    
    # Normalização Intensiva da Matriz Laplaciana
    k_mean_mech = max(1.0, mean(degree(g_mech)))
    k_mean_elec = max(1.0, mean(degree(g_elec)))
    
    L1 = sparse(Matrix{Float64}(laplacian_matrix(g_mech)) ./ k_mean_mech)
    L2 = sparse(Matrix{Float64}(laplacian_matrix(g_elec)) ./ k_mean_elec)

    S_mech = bianconi_entropy(g_mech)
    S_elec = bianconi_entropy(g_elec)
    
    open("$folder/entropia_estrutural_N_$N.txt", "w") do io
        write(io, "N\tBianconi_Mech\tBianconi_Elec\n")
        write(io, "$N\t$S_mech\t$S_elec\n")
    end

    plot_data = Dict{Float64, Matrix{Float64}}()

    for α in alphas
        println("  -> Processando Alpha = $α ...")
        
        resultados_sigma = zeros(length(sigma_inter_range), 5)
        
        # O processamento paralelo distribui os 100 pontos da curva pelos núcleos do computador
        Threads.@threads for i in 1:length(sigma_inter_range)
            σ_inter = sigma_inter_range[i]
            
            # f0 original era 1.2. Aumentado em 20% resulta em f0 = 1.44
            p_phys = (0.1, 0.5, 1.0, 0.5, 0.5, 1, 1.44, 0.8, 0.5, 0.05, σ_inter)
            
            x, v, V_el = solve_piezo_duplex_fractional(N, L1, L2, α, p_phys, t_steps, dt)
            
            # Janela Estacionária: Pegamos os últimos 2 segundos da física (2000 passos)
            idx = t_steps-1999:t_steps
            phases_mech = atan.(v[:, idx], x[:, idx])
            
            R_mech = mean([abs(mean(exp.(im .* phases_mech[:, t]))) for t in 1:size(phases_mech, 2)])
            P_avg = mean(V_el[:, idx].^2)
            H, MI = calculate_info_metrics(phases_mech)
            
            resultados_sigma[i, :] = [σ_inter, R_mech, P_avg, H, MI]
        end
        
        header = "Sigma_Inter\tKuramoto_R\tPower_Avg\tShannon_H\tMutual_Info"
        filename = "$folder/dados_alpha_$(α)_N_$(N).txt"
        writedlm(filename, vcat(permutedims(split(header, "\t")), resultados_sigma), '\t')
        
        plot_data[α] = resultados_sigma
    end

    # --- 5. Geração do Gráfico Comparativo (Dashboard) ---
    println("  -> Gerando gráficos para N = $N ...")
    
    default(fontfamily="Arial", grid=true, gridalpha=0.3, frame=:box, lw=2.5, dpi=300)
    colors = [:blue, :orange, :purple]
    
    p1 = plot(title="Mechanical Synchrony (R)", ylabel="R")
    p2 = plot(title="Harvested Power (P)", ylabel="P")
    p3 = plot(title="Shannon Entropy (H)", ylabel="H", xlabel="Piezo Coupling (σ_inter)")
    p4 = plot(title="Mutual Information (MI)", ylabel="MI", xlabel="Piezo Coupling (σ_inter)")
    
    for (i, α) in enumerate(alphas)
        d = plot_data[α]
        plot!(p1, d[:, 1], d[:, 2], label="α=$α", color=colors[i])
        plot!(p2, d[:, 1], d[:, 3], label="α=$α", color=colors[i])
        plot!(p3, d[:, 1], d[:, 4], label="α=$α", color=colors[i], ls=:solid)
        plot!(p4, d[:, 1], d[:, 5], label="α=$α", color=colors[i], ls=:solid)
    end
    
    painel = plot(p1, p2, p3, p4, layout=(2,2), size=(1000, 800), plot_title="Duplex Energy Harvesting Analysis (N=$N)", margin=5Plots.mm)
    savefig(painel, "$folder/comparativo_metricas_N_$N.png")
end

println("\nTodas as simulações concluídas. Verifique as pastas geradas!")