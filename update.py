import os
import subprocess
from datetime import datetime

def generate_directory_tree(root_dir, level=0, prefix=''):
    """递归生成目录树"""
    extra_extensions = ['.ipynb']  # 在这里添加额外需要包含的文件类型
    
    content = [] if level > 0 else ["# Repository Directory Structure\n\n"]
    has_valid_files = False
    
    items = os.listdir(root_dir)
    files = [f for f in items if os.path.isfile(os.path.join(root_dir, f)) 
             and (f.endswith('.md') or any(f.endswith(ext) for ext in extra_extensions))
             and not (f == 'README.md' and os.path.samefile(root_dir, current_dir))]
    
    if files:
        has_valid_files = True
        for file in sorted(files):
            rel_path = os.path.join(prefix, file) if prefix else file
            content.append(f"{' ' * (level * 4)}- [{os.path.splitext(file)[0]}]({rel_path})\n")
    
    dirs = [d for d in items if os.path.isdir(os.path.join(root_dir, d)) and d != '.git']
    for dir_name in sorted(dirs):
        dir_path = os.path.join(root_dir, dir_name)
        new_prefix = os.path.join(prefix, dir_name) if prefix else dir_name
        
        # 递归处理子目录
        sub_content, sub_has_files = generate_directory_tree(dir_path, level + 1, new_prefix)
        
        # 只有当子目录包含有效文件时，才添加目录名和子目录内容
        if sub_has_files:
            has_valid_files = True
            content.append(f"{' ' * (level * 4)}- {dir_name}/\n")
            content.extend(sub_content)
    
    return content, has_valid_files

def git_commit():
    # 获取当前时间
    current_time = datetime.now().strftime("%Y-%m-%d %H:%M")
    commit_message = f"update in {current_time}"
    try:
        subprocess.run(['git', 'add', '.'], check=True)
        subprocess.run(['git', 'commit', '-m', commit_message], check=True)
        
        print(f"Successfully committed with message: {commit_message}")
    except subprocess.CalledProcessError as e:
        print(f"Error during git operations: {e}")

if __name__ == "__main__":
    current_dir = os.path.dirname(os.path.abspath(__file__))
    content, _ = generate_directory_tree(current_dir)
    
    with open(os.path.join(current_dir, 'README.md'), 'w', encoding='utf-8') as f:
        f.writelines(content)
    
    git_commit()