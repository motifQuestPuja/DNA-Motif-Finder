import streamlit as st
from collections import Counter

# Set up the page
st.set_page_config(page_title="DNA Motif Finder", layout="wide")

# Sidebar - Tool Options
st.sidebar.title("🧬 DNA Motif Finder")
st.sidebar.markdown("### 🔧 Tool Options")

# 1. Input field
motif_input = st.sidebar.text_area("Enter DNA Sequences (FASTA or plain)", height=200)

# 2. Motif length slider (4 to 20)
motif_length = st.sidebar.slider("Select Motif Length", min_value=4, max_value=20, value=6)

# 3. Algorithm selector
algorithm = st.sidebar.selectbox("Select Algorithm", ["Greedy Search"])  # Only one implemented

# 4. Run button
run_button = st.sidebar.button("Run Motif Finder")

st.sidebar.markdown("---")
st.sidebar.markdown("### 👁️ View Options")

# 5. View toggles
view_alignment = st.sidebar.checkbox("Show Sequence Alignment", value=True)
view_consensus = st.sidebar.checkbox("Show Consensus Motif", value=True)
view_pwm = st.sidebar.checkbox("Show Position Weight Matrix", value=True)
view_logo = st.sidebar.checkbox("Show Sequence Logo", value=True)

# Title
st.title("🔍 DNA Motif Finder")

# HOME PAGE tabs
home_tabs = st.tabs(["HOME PAGE"])

with home_tabs[0]:
    about_tab, team_tab = st.tabs(["About Page", "Team Page"])

    with about_tab:
        st.header("About DNA Motif Finder")
        st.write("""
        **DNA Motif Finder** is a web-based application designed to identify conserved motifs within a set of DNA sequences.

        ### Objectives:
        - Identify overrepresented motifs using probabilistic and deterministic approaches.
        - Provide biological insight into potential transcription factor binding sites.
        - Enable interpretation through visual and numerical outputs.

        ### Features:
        - Input multiple DNA sequences in FASTA/plain format.
        - Select from common motif search algorithms.
        - View alignment, consensus motifs, PWM, and sequence logos.
        """)

    with team_tab:
        st.header("Team & Acknowledgement")

        st.subheader("Contributor")
        st.markdown("""
        **Pooja Swamy**  
        Masters in Bioinformatics  
        **University:** DES Pune University School of Science and Mathematics  
        **Skills:** Python, Bioinformatics tools and databases  
        **Project:** Developed DNA Motif Finder as part of academic research project.
        """)

        st.subheader("Guide")
        st.markdown("""
        **Dr. Kushagra Kashyap**  
        Assistant Professor, DES Pune University  
        🔗 [LinkedIn: Dr. Kushagra Kashyap](https://www.linkedin.com)  
        ✉️ Email: 3522411035@despu.edu.in  
        💻 GitHub: *(GitHub link not provided)*
        """)

        st.subheader("Acknowledgement")
        st.write("""
        I would like to express my heartfelt gratitude to **Dr. Kushagra Kashyap** and **Dr. Poonam Deshpande Ma'am** for their invaluable guidance and support throughout the development of the DNA Motif Finder web application.

        This project was undertaken as a **bioinformatics project** and has been a great learning experience in integrating **computational knowledge** with **biological analysis**.

        I am also thankful for the knowledge gained through **online databases, lectures**, and resources that enabled the development of this application.

        **Thank you.**
        """)

# --- Algorithm Functions ---

def count_matrix(motifs):
    k = len(motifs[0])
    count = {'A':[0]*k, 'C':[0]*k, 'G':[0]*k, 'T':[0]*k}
    for motif in motifs:
        for i, nuc in enumerate(motif):
            count[nuc][i] += 1
    return count

def profile_matrix(motifs):
    count = count_matrix(motifs)
    t = len(motifs)
    profile = {nuc: [count[nuc][i] / t for i in range(len(motifs[0]))] for nuc in "ACGT"}
    return profile

def consensus_motif(motifs):
    count = count_matrix(motifs)
    consensus = ""
    for i in range(len(motifs[0])):
        max_freq = 0
        max_nuc = ""
        for nuc in "ACGT":
            if count[nuc][i] > max_freq:
                max_freq = count[nuc][i]
                max_nuc = nuc
        consensus += max_nuc
    return consensus

def score_motifs(motifs):
    consensus = consensus_motif(motifs)
    return sum(sum(1 for i in range(len(motif)) if motif[i] != consensus[i]) for motif in motifs)

def most_probable_kmer(seq, k, profile):
    max_prob = -1
    best_kmer = seq[0:k]
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i+k]
        prob = 1
        for j, nuc in enumerate(kmer):
            prob *= profile[nuc][j]
        if prob > max_prob:
            max_prob = prob
            best_kmer = kmer
    return best_kmer

def greedy_motif_search(dna_seqs, k):
    best_motifs = [seq[0:k] for seq in dna_seqs]
    first_seq = dna_seqs[0]
    
    for i in range(len(first_seq) - k + 1):
        motifs = [first_seq[i:i+k]]
        for j in range(1, len(dna_seqs)):
            profile = profile_matrix(motifs)
            next_motif = most_probable_kmer(dna_seqs[j], k, profile)
            motifs.append(next_motif)
        if score_motifs(motifs) < score_motifs(best_motifs):
            best_motifs = motifs
            
    return best_motifs, consensus_motif(best_motifs), profile_matrix(best_motifs)

# --- Run the motif finder if button is clicked ---
if run_button:
    sequences = [line.strip() for line in motif_input.strip().splitlines() if line and not line.startswith(">")]

    if len(sequences) < 2:
        st.warning("⚠️ Please enter at least two DNA sequences.")
    elif any(len(seq) < motif_length for seq in sequences):
        st.warning("⚠️ All sequences must be longer than the selected motif length.")
    else:
        st.subheader(f"🔬 Running {algorithm} for motif length {motif_length}...")
        motifs, consensus, pwm = greedy_motif_search(sequences, motif_length)
        st.success("✅ Motif search completed!")

        if view_alignment:
            st.subheader("🧬 Motif Alignment")
            for i, motif in enumerate(motifs):
                st.text(f"Seq {i+1}: {motif}")

        if view_consensus:
            st.subheader("🎯 Consensus Motif")
            st.code(consensus)

        if view_pwm:
            st.subheader("📊 Position Weight Matrix")
            st.write(pwm)

        if view_logo:
            st.subheader("🖼️ Sequence Logo")
            st.info("Sequence logo visualization coming soon!")
