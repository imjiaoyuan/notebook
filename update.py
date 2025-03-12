import os
import subprocess
from datetime import datetime

# 配置参数
CONFIG = {
    # 需要包含的文件扩展名
    'include_extensions': [
        '.md',
        '.ipynb',
        '.r',
        '.py',
        '.sh'
    ],
    # 需要排除的文件
    'exclude_files': [
        'update.py',
        'README.md'
    ],
    # 需要排除的目录
    'exclude_dirs': [
        '.git',
        'Dotfiles',
        'Dissertation',
        'ScriptForge',
        'CV'
    ]
}

def normalize_path(path):
    """统一路径分隔符为正斜杠，确保跨平台兼容性"""
    return path.replace(os.path.sep, '/')

def generate_directory_tree(root_dir, level=0):
    """递归生成目录树"""
    content = [] if level > 0 else ["# Repository Directory Structure\n\n"]
    base_dir = os.path.dirname(root_dir)  # 获取父目录
    has_valid_files = False
    
    try:
        items = os.listdir(root_dir)
    except PermissionError:
        return content, False
    
    # 处理文件
    files = [f for f in items if os.path.isfile(os.path.join(root_dir, f)) 
            and (any(f.endswith(ext) for ext in CONFIG['include_extensions']))
            and f not in CONFIG['exclude_files']
            and not (f == 'README.md' and os.path.samefile(root_dir, base_dir))]
    
    if files:
        has_valid_files = True
        for file in sorted(files):
            # 构建相对于仓库根目录的路径
            rel_dir = os.path.relpath(root_dir, base_dir)
            file_path = normalize_path(os.path.join(rel_dir, file))
            content.append(f"- [{file}]({file_path})\n")
    
    # 处理目录
    dirs = [d for d in items if os.path.isdir(os.path.join(root_dir, d)) 
           and d not in CONFIG['exclude_dirs']]
    
    for dir_name in sorted(dirs):
        dir_path = os.path.join(root_dir, dir_name)
        sub_content, sub_has_files = generate_directory_tree(dir_path, level + 1)
        
        if sub_has_files:
            has_valid_files = True
            # 只对第一级目录加粗
            if level == 0:
                content.append(f"- **{dir_name}/**\n")
            else:
                content.append(f"- {dir_name}/\n")
            content.extend([f"  {line}" for line in sub_content])  # 添加缩进
    
    return content, has_valid_files

def git_commit():
    """提交更改到git仓库"""
    try:
        # 获取当前时间
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        
        # 添加所有更改
        subprocess.run(["git", "add", "."], check=True)
        
        # 提交更改
        commit_message = f"Update directory structure: {current_time}"
        subprocess.run(["git", "commit", "-m", commit_message], check=True)
        
        # 可选：推送到远程仓库
        # subprocess.run(["git", "push"], check=True)
        
        print(f"Successfully committed changes with message: {commit_message}")
    except subprocess.CalledProcessError as e:
        print(f"Error during git operations: {str(e)}")
        print("Make sure you have git installed and repository initialized")

if __name__ == "__main__":
    current_dir = os.path.dirname(os.path.abspath(__file__))
    content, _ = generate_directory_tree(current_dir)
    
    readme_path = os.path.join(current_dir, 'README.md')
    with open(readme_path, 'w', encoding='utf-8') as f:
        f.writelines(content)
    
    print(f"Directory structure has been written to {readme_path}")
    
    git_commit()