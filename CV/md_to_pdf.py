from markdown import markdown
import jinja2
from weasyprint import HTML, CSS
from pathlib import Path
import re

def read_markdown_file(file_path):
    """读取Markdown文件内容"""
    with open(file_path, 'r', encoding='utf-8') as f:
        return f.read()

def convert_markdown_to_html(md_content):
    """将Markdown转换为HTML，并进行自定义处理"""
    # 移除HTML注释
    md_content = re.sub(r'<!--.*?-->', '', md_content, flags=re.DOTALL)
    
    # 转换Markdown为HTML
    html_content = markdown(md_content, extensions=['extra'])
    
    # 处理自定义的hr样式，移除阴影效果
    html_content = html_content.replace(
        '<hr>',
        '<hr style="border: none; border-top: 1px solid #000000; margin: 6px 0;">'
    )
    
    # 增强加粗效果，移除阴影
    bold_style = (
        'style="'
        'color: #000000; '
        'font-weight: 900; '
        'letter-spacing: 0.02em; '
        'display: inline-block; '
        '-webkit-font-smoothing: antialiased; '
        '-moz-osx-font-smoothing: grayscale; '
        'text-rendering: optimizeLegibility; '
        '"'
    )
    html_content = html_content.replace('<strong>', f'<strong {bold_style}>')
    
    return html_content

def get_html_template():
    """获取HTML模板"""
    return """
    <!DOCTYPE html>
    <html>
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0, user-scalable=yes">
        <title>Resume</title>
        <style>
            @page {
                size: A4;
                margin: 0.4cm 1cm 1cm 1cm; /*上 右 下 左*/
            }
            
            @media print {
                body {
                    margin: 0;
                    padding: 0;
                }
            }
            
            body {
                font-family: "Noto Sans CJK SC", "Microsoft YaHei", sans-serif;
                font-size: 14px;
                line-height: 1.5;
                color: #333333;
                word-wrap: break-word;
                max-width: 21cm;
                margin: 0 auto;
                padding: 0.4em;
                background: #fff;
                -webkit-font-smoothing: antialiased;  /* 增加字体清晰度 */
                -moz-osx-font-smoothing: grayscale;
                text-rendering: optimizeLegibility;
            }
            
            h3 {
                color: #000000;
                font-weight: 900;  /* 使用最大字重 */
                padding-bottom: 0.1em;
                font-size: 1.2em;
                border-bottom: none;
                margin-top: 10px;
                margin-bottom: 6px;
                line-height: 1.2;
            }
            
            hr {
                border: none;
                border-top: 1px solid #000000;
                margin: 6px 0;
                padding: 0;
                background: none;  /* 移除背景色 */
                height: 0;  /* 确保只有边框线 */
            }
            
            h3 + ul {
                margin-top: 6px;
            }
            
            ul {
                margin-top: 0;
                margin-bottom: 8px;
                padding-left: 2em;
                line-height: 1.5;
            }
            
            li {
                margin: 0;
                padding: 0;
                line-height: 1.5;
            }
            
            /* 调整列表间距 */
            li ul {
                margin-top: 0;  /* 子列表顶部间距 */
                margin-bottom: 0;  /* 子列表底部间距 */
            }
            
            li li {
                margin-top: 0;  /* 子列表项间距 */
            }
            
            /* 调整相邻列表项间距 */
            li + li {
                margin-top: 0;  /* 移除相邻列表项的额外间距 */
            }
            
            /* 加粗文本样式 */
            strong, b {
                color: #000000;
                font-weight: 900;
                letter-spacing: 0.02em;
                display: inline-block;
                -webkit-font-smoothing: antialiased;
                -moz-osx-font-smoothing: grayscale;
                text-rendering: optimizeLegibility;
            }
            
            p {
                margin-top: 0;
                margin-bottom: 8px;  /* 减小段落间距 */
            }
            
            a {
                color: #0969da;
                text-decoration: none;
            }
            
            a:hover {
                text-decoration: underline;
            }
            
            img {
                max-width: 80px;
                float: right;
                margin-right: 25px;
            }
            
            code {
                padding: .2em .4em;
                margin: 0;
                font-size: 85%;
                background-color: rgba(175,184,193,0.2);
                border-radius: 6px;
                font-family: Consolas, monospace;
            }
            
            pre {
                padding: 16px;
                overflow: auto;
                font-size: 85%;
                line-height: 1.45;
                background-color: #f6f8fa;
                border-radius: 6px;
            }
            
            pre code {
                padding: 0;
                margin: 0;
                background-color: transparent;
            }
        </style>
    </head>
    <body>
        {{ content }}
    </body>
    </html>
    """

def generate_pdf(md_file_path, output_path):
    """生成PDF文件"""
    # 读取Markdown内容
    md_content = read_markdown_file(md_file_path)
    
    # 转换为HTML
    html_content = convert_markdown_to_html(md_content)
    
    # 渲染HTML模板
    template = jinja2.Template(get_html_template())
    html_output = template.render(content=html_content)
    
    # 生成PDF
    html = HTML(string=html_output)
    css = CSS(string='''
        @page {
            size: A4;
            margin: 0.8cm 1cm 1.2cm 1cm;  /* 确保与模板中的边距一致 */
        }
    ''')
    
    # 确保输出目录存在
    output_dir = Path(output_path).parent
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # 生成PDF
    html.write_pdf(output_path, stylesheets=[css])
    print(f"PDF已生成: {output_path}")

if __name__ == '__main__':
    # 设置输入输出路径
    current_dir = Path(__file__).parent
    input_file = current_dir / 'CV.md'
    output_file = current_dir / 'CV.pdf'
    
    # 生成PDF
    generate_pdf(input_file, output_file) 