from markdown import markdown
import jinja2
from weasyprint import HTML, CSS
from pathlib import Path
import re

def read_markdown_file(file_path):
    with open(file_path, 'r', encoding='utf-8') as f:
        return f.read()

def convert_markdown_to_html(md_content):
    md_content = re.sub(r'<!--.*?-->', '', md_content, flags=re.DOTALL)
    html_content = markdown(md_content, extensions=['extra', 'tables', 'sane_lists'])
    
    html_content = html_content.replace('<hr>', '<hr style="border: none; border-top: 1px solid #000000; margin: 6px 0;">')
    
    bold_style = 'style="color: #000000; font-weight: 900; letter-spacing: 0.02em;"'
    return html_content.replace('<strong>', f'<strong {bold_style}>')

def get_html_template():
    return """<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <title>Resume</title>
    <style>
        @page {
            size: A4;
            margin: 0.5cm 1cm;
        }
        
        body {
            font-family: "Noto Sans CJK SC", "Microsoft YaHei", sans-serif;
            font-size: 12pt;
            line-height: 1.5;
            color: #333;
            max-width: 21cm;
            margin: 0 auto;
            padding: 0.5cm;
        }
        
        h3 {
            color: #000;
            font-weight: bold;
            font-size: 1.2em;
            margin: 12px 0 6px 0;
            padding-bottom: 2px;
        }
        
        hr {
            border: none;
            border-top: 1px solid #000;
            margin: 6px 0;
        }
        
        ul, ol {
            margin: 4px 0;
            padding-left: 1.5em;
        }
        
        li {
            margin: 4px 0;
            padding: 0;
        }
        
        ul ul, ol ol {
            margin: 2px 0;
            padding-left: 1.5em;
        }
        
        strong, b {
            color: #000;
            font-weight: bold;
        }
        
        p {
            margin: 6px 0;
        }
    </style>
</head>
<body>
    {{ content }}
</body>
</html>"""

def generate_pdf(md_file_path, output_path):
    html_content = convert_markdown_to_html(read_markdown_file(md_file_path))
    template = jinja2.Template(get_html_template())
    html_output = template.render(content=html_content)
    
    base_url = f"file://{str(Path(md_file_path).parent.absolute())}/"
    html = HTML(string=html_output, base_url=base_url)
    
    output_dir = Path(output_path).parent
    output_dir.mkdir(parents=True, exist_ok=True)
    
    html.write_pdf(output_path)
    print(f"PDF已生成: {output_path}")

if __name__ == '__main__':
    current_dir = Path(__file__).parent
    input_file = current_dir / '../CV/CV.md'
    output_file = current_dir / '../CV/CV.pdf'
    generate_pdf(input_file, output_file)