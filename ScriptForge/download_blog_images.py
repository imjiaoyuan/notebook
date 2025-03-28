import os
import re
import requests
from urllib.parse import urlparse
from pathlib import Path
import hashlib
import concurrent.futures
from tqdm import tqdm
import time
from datetime import datetime

class ImageDownloader:
    def __init__(self, posts_dir, image_dir, max_retries=3, max_workers=5):
        self.posts_dir = posts_dir
        self.image_dir = image_dir
        self.max_retries = max_retries
        self.max_workers = max_workers
        self.failed_urls_file = "failed_downloads.txt"
        self.download_log_file = "download_log.txt"
        self.pattern = r'https://images\.yuanj\.top/blog/[^\s\)\"\']+\.(?:png|jpg|jpeg|gif)'
        
        # 设置代理
        self.proxies = {
            'http': 'http://127.0.0.1:7897',
            'https': 'http://127.0.0.1:7897'
        }
        
        # 确保目录存在
        os.makedirs(image_dir, exist_ok=True)
        
        # 初始化统计信息
        self.stats = {
            'total_images': 0,
            'downloaded': 0,
            'skipped': 0,
            'failed': 0,
            'start_time': None,
            'end_time': None
        }
        
        # 初始化集合
        self.downloaded_hashes = set()
        self.all_image_urls = []
        self.failed_urls = []

    def log_message(self, message):
        """记录日志"""
        timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        log_message = f"[{timestamp}] {message}"
        print(log_message)
        with open(self.download_log_file, 'a', encoding='utf-8') as f:
            f.write(log_message + '\n')

    def download_image(self, url):
        """下载单个图片"""
        for attempt in range(self.max_retries):
            try:
                response = requests.get(url, timeout=10, proxies=self.proxies)
                if response.status_code == 200:
                    filename = os.path.basename(urlparse(url).path)
                    save_path = os.path.join(self.image_dir, filename)
                    
                    content_hash = hashlib.md5(response.content).hexdigest()
                    
                    if content_hash in self.downloaded_hashes:
                        return {"status": "skipped", "filename": filename}
                    
                    with open(save_path, 'wb') as f:
                        f.write(response.content)
                    
                    self.downloaded_hashes.add(content_hash)
                    return {"status": "success", "filename": filename}
                else:
                    if attempt == self.max_retries - 1:
                        return {"status": "failed", "filename": url, 
                               "error": f"HTTP {response.status_code}"}
            except Exception as e:
                if attempt == self.max_retries - 1:
                    return {"status": "failed", "filename": url, "error": str(e)}
                self.log_message(f"重试 {attempt + 1}/{self.max_retries}: {url}")
                time.sleep(1)  # 重试前等待1秒
        return {"status": "failed", "filename": url, "error": "达到最大重试次数"}

    def collect_image_urls(self):
        """收集所有markdown文件中的图片URL"""
        self.log_message(f"开始扫描目录: {self.posts_dir}")
        
        for root, _, files in os.walk(self.posts_dir):
            for file in files:
                if file.endswith('.md'):
                    file_path = os.path.join(root, file)
                    self.log_message(f"处理文件: {file_path}")
                    
                    try:
                        with open(file_path, 'r', encoding='utf-8') as f:
                            content = f.read()
                        
                        image_urls = re.findall(self.pattern, content)
                        self.all_image_urls.extend(image_urls)
                        self.stats['total_images'] += len(image_urls)
                        
                    except Exception as e:
                        self.log_message(f"处理文件失败 {file}: {str(e)}")

    def download_all_images(self):
        """并行下载所有图片"""
        self.stats['start_time'] = datetime.now()
        self.log_message("开始下载图片...")
        
        with concurrent.futures.ThreadPoolExecutor(max_workers=self.max_workers) as executor:
            future_to_url = {
                executor.submit(self.download_image, url): url 
                for url in self.all_image_urls
            }
            
            with tqdm(total=len(self.all_image_urls), desc="下载进度") as pbar:
                for future in concurrent.futures.as_completed(future_to_url):
                    url = future_to_url[future]
                    result = future.result()
                    
                    if result["status"] == "success":
                        self.stats['downloaded'] += 1
                        self.log_message(f"✅ 成功下载: {result['filename']}")
                    elif result["status"] == "skipped":
                        self.stats['skipped'] += 1
                        self.log_message(f"⏭️ 跳过重复图片: {result['filename']}")
                    else:
                        self.stats['failed'] += 1
                        self.failed_urls.append(url)
                        self.log_message(f"❌ 下载失败 {url}: {result['error']}")
                    pbar.update(1)
        
        self.stats['end_time'] = datetime.now()

    def save_failed_urls(self):
        """保存失败的URL到文件"""
        if self.failed_urls:
            with open(self.failed_urls_file, 'w', encoding='utf-8') as f:
                for url in self.failed_urls:
                    f.write(f"{url}\n")
            self.log_message(f"失败的URL已保存到 {self.failed_urls_file}")

    def print_statistics(self):
        """打印统计信息"""
        duration = self.stats['end_time'] - self.stats['start_time']
        
        self.log_message("\n📊 下载统计:")
        self.log_message(f"开始时间: {self.stats['start_time'].strftime('%Y-%m-%d %H:%M:%S')}")
        self.log_message(f"结束时间: {self.stats['end_time'].strftime('%Y-%m-%d %H:%M:%S')}")
        self.log_message(f"总耗时: {duration}")
        self.log_message(f"总共发现图片: {self.stats['total_images']}")
        self.log_message(f"成功下载图片: {self.stats['downloaded']}")
        self.log_message(f"跳过重复图片: {self.stats['skipped']}")
        self.log_message(f"失败数量: {self.stats['failed']}")

    def run(self):
        """运行下载流程"""
        self.collect_image_urls()
        if not self.all_image_urls:
            self.log_message("没有找到需要下载的图片")
            return
        
        self.download_all_images()
        self.save_failed_urls()
        self.print_statistics()

def main():
    # 设置路径
    posts_dir = r"D:\Code\imjiaoyuan.github.io\content\posts"
    image_dir = r"D:\Code\imjiaoyuan.github.io\content\downloaded_images"
    
    # 检查目录是否存在
    if not os.path.exists(posts_dir):
        print(f"❌ 错误: 文章目录不存在: {posts_dir}")
        exit(1)
    
    # 创建下载器实例并运行
    downloader = ImageDownloader(
        posts_dir=posts_dir,
        image_dir=image_dir,
        max_retries=3,
        max_workers=5
    )
    downloader.run()

if __name__ == "__main__":
    main() 