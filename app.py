#Imports
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
import gradio as gr
import time
from google import genai
from google.genai import types


def blast_search(file_path, queries, alg, api):
    client=genai.Client(api_key=api)
    
    text=""
    #Takes the contents of the file into a string
    fasta_string=open(file_path.name).read()
    
    #Splits the file into records
    records=fasta_string.split(">")
    
    #Removes the blank string
    records.remove(records[0])
    
    #Re-adds the ">" sign in each record
    for record in records:
        records[records.index(record)] = ">"+record
    
    queries = queries.split(",")
    
    #Turns the strings to integers
    for query in queries:
        queries[queries.index(query)]=int(query)
    
    alg=alg.lower()
    if alg=="blastn":
        mb=None
    elif alg=="megablast":
        mb=True
        
    #Searches using qblast
    for query in queries:
        start_time=time.perf_counter()
        text+=("\n\n"*3)
        text+= f"**** Searching Record Number {query} ****"
        yield text
        result_handle=NCBIWWW.qblast("blastn","nt",records[query-1], megablast=mb)
        
        #Reads these results into the record
        blast_record=NCBIXML.read(result_handle)
    
        #Sets the e value threshold (This program uses 0.01, default for BLAST is 0.05), so this is MORE selective
        e_value_thresh=0.01
        
        #Checks if there are alignments to print
        #Shows alignments

        yield text+"\n\n"+"Generating Gemini Report on Organism"
        top_alignment_report=data=client.models.generate_content(
        model='gemini-2.5-pro',
        config=types.GenerateContentConfig(
        temperature=0.0,
        top_p=1.0),
        contents=[f"Please generate a short but descriptive summary of the following organism - {blast_record.alignments[0].title}"]).text
        text+=f"\n\nRecord's Top Alignment Report - {top_alignment_report}"
        yield text

        
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                if hsp.expect < e_value_thresh:
                    text+="\n\n**Alignment**"+f"\n\nSequence: **{alignment.title.split('| ')[1]}**"+f"\n\nSequence Identity: {alignment.title.split('| ')[0]}"+f"\n\nLength: {alignment.length}"+f"\n\nE Value: {hsp.expect}"+f"\n\nScore: {hsp.score}"
                    yield text
                    if hsp.query==hsp.sbjct:
                        text+= "\n\n**Exact Match!**"
                        end_time=time.perf_counter()
                        execution_time=abs(end_time-start_time)
                        text+= f"\n\nExecution time: {round(execution_time)} seconds"
                        text+= "\n\n"
                        yield text
                    else:
                        text+= "\n\n"+str(hsp.query)
                        text+= "\n\n"+str(hsp.match)
                        text+= "\n\n"+str(hsp.sbjct)
                        end_time=time.perf_counter()
                        execution_time=abs(end_time-start_time)
                        text+= f"\n\nExecution time: {round(execution_time)} seconds"
                        yield text+ "\n\n"
                        
    #Finishes results
    text+= "\n\n No other significant similarities found"
    yield text

with gr.Blocks(theme=gr.themes.Soft(), title="Biopython BLAST Search") as demo:
    gr.Markdown("# 🧬 Biopython BLAST Search Engine")
    gr.Markdown("Provide your sequence(s) in FASTA format below.")

    with gr.Row():
        with gr.Column(scale=1):
            # Use Tabs to create separate input areas
            with gr.Tabs():
                with gr.TabItem("Upload File"):
                    fasta_file_input = gr.File(label="Upload any file with FASTA sequences (.fa, .txt, etc.)")
            query_input = gr.Textbox(label="Record numbers to search", placeholder="e.g., 1, 3, 5")
            algorithm_input = gr.Radio(choices=["blastn", "megablast"], label="Search Algorithm", value="blastn")
            api_key=gr.Textbox(label='Gemini API Key for Reports', type="password")
            submit_button = gr.Button("Run BLAST Search", variant="primary")
        with gr.Column(scale=2):
            results_output = gr.Markdown("Your results will appear here...")

    # The click event now passes both input components to the function
    submit_button.click(
        fn=blast_search,
        inputs=[fasta_file_input, query_input, algorithm_input, api_key],
        outputs=results_output
    )

# Launch the Gradio app
demo.launch()